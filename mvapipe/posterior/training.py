import csv
import json
import random
import statistics
import zipfile
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from enum import Enum
from pathlib import Path
from typing import TypedDict

import numpy as np

from mvapipe.frequency import FrequencyMode
from mvapipe.gene_validity import UNKNOWN_GENE, GeneValidity, Inheritance, Validity
from mvapipe.ontology import Ontology, PhenotypeModel
from mvapipe.posterior.folds import read_gene_folds
from mvapipe.posterior.phenotype import (
    LeaveSourceOut,
    PhenotypeTable,
    SimilarityIndex,
    find_pubmed_ids,
    read_annotation_sources,
)
from mvapipe.posterior.priors import Genotype
from mvapipe.references import References
from mvapipe.variant import Variant
from mvapipe.variant_table import load_variants
from mvapipe.vcf_steps import CANONICAL_CHROMOSOMES

SEED = 20261018
NEGATIVE_GENES_PER_CASE = 200
AUTOSOMES = tuple(CANONICAL_CHROMOSOMES[:22])
ALLELIC_STATES = {
    "GENO:0000135": "het",
    "GENO:0000136": "hom",
    "GENO:0000134": "hemi",
}
MAXIMUM_ALLELE_LENGTH = 51


class Context(Enum):
    """Genotype contexts fitted with separate variant models."""

    DOMINANT = "dominant"
    RECESSIVE_HETEROZYGOUS = "recessive_heterozygous"
    RECESSIVE_HOMOZYGOUS = "recessive_homozygous"


CONTEXT_OF_GENOTYPE = {
    Genotype.DOMINANT: Context.DOMINANT,
    Genotype.COMPOUND: Context.RECESSIVE_HETEROZYGOUS,
    Genotype.SINGLE_RECESSIVE: Context.RECESSIVE_HETEROZYGOUS,
    Genotype.HOMOZYGOUS: Context.RECESSIVE_HOMOZYGOUS,
}


@dataclass(frozen=True)
class SmallAllele:
    # Key format: chromosome-position-reference-alternate, without the 'chr' prefix.
    key: str
    state: str
    gene: str


@dataclass(frozen=True)
class PublishedCase:
    identifier: str
    terms: tuple[str, ...]
    small_alleles: tuple[SmallAllele, ...]
    # Causal alleles without a GRCh38 small-variant record.
    other_allele_count: int
    pubmed_ids: frozenset[str]
    other_genes: frozenset[str] = frozenset()
    has_unknown_causal_gene: bool = False


# Only the Phenopacket fields used by this parser are declared below.
class OntologyClass(TypedDict, total=False):
    id: str


class VcfRecord(TypedDict, total=False):
    genomeAssembly: str
    chrom: int | str
    pos: int | str
    ref: str
    alt: str


class GeneDescriptor(TypedDict, total=False):
    symbol: str


class VariationDescriptor(TypedDict, total=False):
    vcfRecord: VcfRecord
    allelicState: OntologyClass
    geneContext: GeneDescriptor


class VariantInterpretation(TypedDict, total=False):
    variationDescriptor: VariationDescriptor


class GenomicInterpretation(TypedDict, total=False):
    interpretationStatus: str
    variantInterpretation: VariantInterpretation


class Diagnosis(TypedDict, total=False):
    genomicInterpretations: list[GenomicInterpretation]


class Interpretation(TypedDict, total=False):
    diagnosis: Diagnosis


class PhenotypicFeature(TypedDict, total=False):
    type: OntologyClass
    excluded: bool


class ExternalReference(TypedDict, total=False):
    id: str


class MetaData(TypedDict, total=False):
    externalReferences: list[ExternalReference]


class Phenopacket(TypedDict, total=False):
    id: str
    metaData: MetaData
    interpretations: list[Interpretation]
    phenotypicFeatures: list[PhenotypicFeature]


def _parse_small_allele(descriptor: VariationDescriptor) -> SmallAllele | None:
    record = descriptor.get("vcfRecord")
    state = ALLELIC_STATES.get(descriptor.get("allelicState", {}).get("id", ""))
    if record is None or state is None:
        return None
    if record.get("genomeAssembly", "").lower() not in ("hg38", "grch38"):
        return None
    chromosome = str(record["chrom"]).removeprefix("chr")
    reference, alternate = record["ref"].upper(), record["alt"].upper()
    is_small = max(len(reference), len(alternate)) <= MAXIMUM_ALLELE_LENGTH
    is_sequence = set(reference + alternate) <= set("ACGT")
    if chromosome not in (*AUTOSOMES, "X") or not is_small or not is_sequence:
        return None
    return SmallAllele(
        key=f"{chromosome}-{int(record['pos'])}-{reference}-{alternate}",
        state=state,
        gene=descriptor.get("geneContext", {}).get("symbol", ""),
    )


def _parse_case(packet: Phenopacket) -> PublishedCase | None:
    small_alleles: list[SmallAllele] = []
    other_allele_count = 0
    other_genes: set[str] = set()
    has_unknown_causal_gene = False
    for interpretation in packet["interpretations"]:
        diagnosis = interpretation.get("diagnosis", {})
        for genomic in diagnosis.get("genomicInterpretations", []):
            if genomic.get("interpretationStatus") != "CAUSATIVE":
                continue
            descriptor = genomic.get("variantInterpretation", {}).get(
                "variationDescriptor"
            )
            allele = None if descriptor is None else _parse_small_allele(descriptor)
            if allele is None:
                other_allele_count += 1
                gene = (
                    descriptor.get("geneContext", {}).get("symbol", "")
                    if descriptor is not None
                    else ""
                )
                if gene:
                    other_genes.add(gene)
                else:
                    has_unknown_causal_gene = True
            else:
                small_alleles.append(allele)
    if not small_alleles and not other_allele_count:
        return None
    features = packet.get("phenotypicFeatures", [])
    return PublishedCase(
        identifier=packet["id"],
        terms=tuple(
            feature["type"]["id"]
            for feature in features
            if not feature.get("excluded", False)
        ),
        small_alleles=tuple(small_alleles),
        other_allele_count=other_allele_count,
        pubmed_ids=frozenset(),
        other_genes=frozenset(other_genes),
        has_unknown_causal_gene=has_unknown_causal_gene,
    )


def read_published_cases(archive_path: Path) -> list[PublishedCase]:
    """Read cases with causal alleles from a Phenopacket Store ZIP archive.

    Read files in filename order and keep the first usable record for each case ID.
    Collect PubMed IDs from every phenopacket sharing that case ID.
    """
    cases: dict[str, PublishedCase] = {}
    pubmed_ids: dict[str, set[str]] = defaultdict(set)
    with zipfile.ZipFile(archive_path) as archive:
        for name in sorted(archive.namelist()):
            if not name.endswith(".json"):
                continue
            packet: Phenopacket = json.loads(archive.read(name))
            if not isinstance(packet, dict) or "id" not in packet:
                continue
            references = packet.get("metaData", {}).get("externalReferences", [])
            pubmed_ids[packet["id"]] |= find_pubmed_ids(
                reference.get("id", "") for reference in references
            )
            if "interpretations" not in packet or packet["id"] in cases:
                continue
            case = _parse_case(packet)
            if case is not None:
                cases[case.identifier] = case
    return [
        replace(case, pubmed_ids=frozenset(pubmed_ids[case.identifier]))
        for case in cases.values()
    ]


def _compute_genomic_position(key: str) -> tuple[int, int, str]:
    chromosome, position, *_ = key.split("-")
    return (*AUTOSOMES, "X").index(chromosome), int(position), key


def write_allele_vcf(keys: Iterable[str], vcf_path: Path) -> None:
    """Write unique alleles in genomic order as heterozygous calls, using keys as IDs."""
    lines = [
        "##fileformat=VCFv4.2",
        *(f"##contig=<ID={chromosome}>" for chromosome in (*AUTOSOMES, "X")),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCASES",
    ]
    for key in sorted(set(keys), key=_compute_genomic_position):
        chromosome, position, reference, alternate = key.split("-")
        lines.append(
            f"{chromosome}\t{position}\t{key}\t{reference}\t{alternate}\t.\tPASS\t.\tGT\t0/1"
        )
    vcf_path.write_text("\n".join(lines) + "\n", encoding="utf8")


@dataclass(frozen=True)
class TrainingCase:
    identifier: str
    gene: str
    genotype: Genotype
    allele_keys: tuple[str, ...]
    # Unique current HPO patient features in phenopacket order.
    terms: tuple[str, ...]
    pubmed_ids: frozenset[str]


def classify_genotype(case: PublishedCase, validity: GeneValidity) -> Genotype | None:
    """Classify a published genotype, returning None for unsupported combinations.

    For a lone heterozygote, require at least limited GenCC evidence to call it dominant.
    If other causal alleles lack small-variant records, classify it as recessive.
    """
    keys = {allele.key for allele in case.small_alleles}
    genes = {allele.gene for allele in case.small_alleles} | case.other_genes
    if case.has_unknown_causal_gene or len(genes) != 1:
        return None
    states = sorted(allele.state for allele in case.small_alleles)
    if case.other_allele_count == 0:
        if states == ["het"]:
            dominant = validity.for_inheritance(Inheritance.AUTOSOMAL_DOMINANT)
            return Genotype.DOMINANT if dominant >= Validity.LIMITED else None
        if states == ["hom"]:
            return Genotype.HOMOZYGOUS
        if states == ["het", "het"] and len(keys) == 2:
            return Genotype.COMPOUND
        return None
    if states == ["het"]:
        return Genotype.SINGLE_RECESSIVE
    return None


def select_training_cases(
    published: Sequence[PublishedCase],
    gene_validity: Mapping[str, GeneValidity],
    ontology: Ontology,
) -> list[TrainingCase]:
    """Select cases with one autosomal gene, a supported genotype and HPO features."""
    cases = []
    for case in published:
        genes = {allele.gene for allele in case.small_alleles} | case.other_genes
        chromosomes = {allele.key.split("-")[0] for allele in case.small_alleles}
        if (
            case.has_unknown_causal_gene
            or "" in genes
            or len(genes) != 1
            or not chromosomes <= set(AUTOSOMES)
        ):
            continue
        gene = genes.pop()
        genotype = classify_genotype(case, gene_validity.get(gene, UNKNOWN_GENE))
        if genotype is None:
            continue
        terms: list[str] = []
        for term in case.terms:
            resolved = ontology.resolve(term)
            if (
                resolved is not None
                and ontology.is_feature(resolved)
                and resolved not in terms
            ):
                terms.append(resolved)
        if not terms:
            continue
        cases.append(
            TrainingCase(
                identifier=case.identifier,
                gene=gene,
                genotype=genotype,
                allele_keys=tuple(sorted({allele.key for allele in case.small_alleles})),
                terms=tuple(terms),
                pubmed_ids=case.pubmed_ids,
            )
        )
    return cases


@dataclass(frozen=True)
class GeneShares:
    """gnomAD coding-variant shares by gene, with a fallback for missing genes."""

    shares: dict[str, float]
    # Use the median share for genes absent from `shares`.
    default_share: float


def read_gene_shares(constraint_path: Path) -> GeneShares:
    """Compute gene shares from gnomAD MANE Select coding-variant counts.

    Sum loss-of-function, missense and synonymous counts. Divide by the total observed
    count, flooring each gene's numerator at 1. Use the median share for missing genes.
    """
    observed: dict[str, float] = {}
    with constraint_path.open(encoding="utf8", newline="") as table:
        for row in csv.DictReader(table, delimiter="\t"):
            if row["mane_select"] != "true":
                continue
            counts = [
                float(row[column])
                for column in ("lof.obs", "mis.obs", "syn.obs")
                if row[column] not in ("", "NA")
            ]
            if counts:
                observed[row["gene"]] = observed.get(row["gene"], 0.0) + sum(counts)
    total = sum(observed.values())
    shares = {gene: max(count, 1.0) / total for gene, count in observed.items()}
    return GeneShares(shares=shares, default_share=statistics.median(shares.values()))


def read_annotated_alleles(
    table_path: Path, *, frequency_mode: FrequencyMode = FrequencyMode.ALL
) -> dict[str, list[Variant]]:
    """Group annotated variants by the allele keys stored in their IDs."""
    alleles: dict[str, list[Variant]] = defaultdict(list)
    for variant in load_variants(table_path, frequency_mode=frequency_mode):
        for identifier in variant.identifier.split(";"):
            alleles[identifier].append(variant)
    return dict(alleles)


def collect_positive_variants(
    cases: Sequence[TrainingCase],
    alleles: Mapping[str, Sequence[Variant]],
    gene_folds: Mapping[str, int],
) -> dict[Context, list[tuple[Variant, int]]]:
    """Collect unique causal alleles per genotype context with their case gene's fold.

    Prefer the annotation for the case's gene; fall back to the allele's first annotation.
    Skip unannotated alleles, including those removed for GRCh38 REF mismatches.
    """
    variants: dict[Context, dict[str, tuple[Variant, int]]] = defaultdict(dict)
    for case in cases:
        context = CONTEXT_OF_GENOTYPE[case.genotype]
        for key in case.allele_keys:
            annotated = alleles.get(key, [])
            if not annotated or key in variants[context]:
                continue
            variant = next((v for v in annotated if v.gene == case.gene), annotated[0])
            variants[context][key] = (variant, gene_folds[case.gene])
    return {context: list(by_key.values()) for context, by_key in variants.items()}


@dataclass(frozen=True)
class TrainingData:
    """Published cases, annotated causal alleles and references for posterior fitting."""

    cases: list[TrainingCase]
    alleles: dict[str, list[Variant]]
    gene_validity: dict[str, GeneValidity]
    similarity: SimilarityIndex
    leave_source_out: LeaveSourceOut
    gene_shares: GeneShares
    gene_folds: Mapping[str, int] | None = None

    @property
    def gene_universe(self) -> list[str]:
        return sorted(self.gene_shares.shares)


def load_training_data(
    references: References,
    gene_validity: dict[str, GeneValidity],
    phenotype_model: PhenotypeModel,
    *,
    allele_table_path: Path | None = None,
    gene_folds_path: Path | None = None,
    frequency_mode: FrequencyMode = FrequencyMode.ALL,
) -> TrainingData:
    """Load training cases, prepared allele annotations and phenotype references."""
    similarity = SimilarityIndex(phenotype_model)
    sources = read_annotation_sources(
        references.phenotype_annotations, phenotype_model.ontology
    )
    published = read_published_cases(references.phenopackets)
    return TrainingData(
        cases=select_training_cases(published, gene_validity, phenotype_model.ontology),
        alleles=read_annotated_alleles(
            allele_table_path or references.training_alleles,
            frequency_mode=frequency_mode,
        ),
        gene_validity=gene_validity,
        similarity=similarity,
        leave_source_out=LeaveSourceOut(similarity, sources),
        gene_shares=read_gene_shares(references.gene_constraint),
        gene_folds=read_gene_folds(gene_folds_path or references.gene_folds),
    )


def build_phenotype_table(
    training: TrainingData, gene_folds: Mapping[str, int], generator: random.Random
) -> PhenotypeTable:
    """Build phenotype training rows for each case's causal gene and sampled genes.

    Draw NEGATIVE_GENES_PER_CASE genes and exclude the causal gene if sampled. Remove
    annotations sourced only from the case's publications. Assign all rows for a case
    to its causal gene's fold.
    """
    gene_index = training.similarity.gene_index
    gene_universe = training.gene_universe
    similarity, is_annotated, is_causal, folds, order = [], [], [], [], []
    for case in training.cases:
        scores, is_gene_annotated = training.leave_source_out.score_genes(
            case.terms, case.pubmed_ids
        )
        sampled = generator.sample(gene_universe, NEGATIVE_GENES_PER_CASE)
        genes = [case.gene] + [gene for gene in sampled if gene != case.gene]
        for position, gene in enumerate(genes):
            row = gene_index.get(gene)
            similarity.append(0.0 if row is None else float(scores[row]))
            is_annotated.append(row is not None and bool(is_gene_annotated[row]))
            is_causal.append(position == 0)
            folds.append(gene_folds[case.gene])
            order.append(position)
    return PhenotypeTable(
        similarity=np.array(similarity),
        is_annotated=np.array(is_annotated, dtype=bool),
        is_causal=np.array(is_causal, dtype=bool),
        fold=np.array(folds),
        order=np.array(order),
    )
