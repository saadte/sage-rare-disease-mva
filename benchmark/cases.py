import functools
import json
import random
import re
from collections import defaultdict
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TypedDict

from mvapipe.gene_validity import UNKNOWN_GENE, GeneValidity, Inheritance, Validity
from mvapipe.ontology import PHENOTYPIC_ABNORMALITY, Ontology, PhenotypeModel
from mvapipe.process import capture
from mvapipe.ranking import parse_review_stars
from mvapipe.variant import Allele

AUTOSOMES = [str(number) for number in range(1, 23)]
PATHOGENIC = frozenset(
    {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
)
LOSS_OF_FUNCTION_CONSEQUENCES = (
    "nonsense",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
)
MINIMUM_REVIEW_STARS = 2
MAXIMUM_ALLELE_LENGTH = 50
CLINVAR_QUERY_FORMAT = (
    "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/GENEINFO\\t%INFO/MC"
    "\\t%INFO/CLNSIG\\t%INFO/CLNREVSTAT\\t%INFO/CLNDISDB\\n"
)
MINIMUM_TERMS_PER_DISEASE = 5
TERMS_PER_CASE = 5
UNRELATED_TERMS_PER_CASE = 2
PARENT_TERM_PROBABILITY = 0.5


class SpikeKind(Enum):
    LOSS_OF_FUNCTION = "lof"
    MISSENSE = "missense"
    OTHER = "other"


class Genotype(Enum):
    COMPOUND_HETEROZYGOUS = "comphet"
    HOMOZYGOUS = "hom"
    HETEROZYGOUS = "het"


class _SpikeRecord(TypedDict):
    chrom: str
    pos: int
    ref: str
    alt: str
    id: str
    gene: str
    kind: str
    diseases: list[str]


class _CaseRecord(TypedDict):
    case: str
    type: str
    gene: str
    variants: list[_SpikeRecord]
    disease: str
    disease_matches_clinvar: bool
    exact: list[str]
    noisy: list[str]


@dataclass(frozen=True)
class SpikeVariant:
    allele: Allele
    clinvar_id: str
    gene: str
    kind: SpikeKind
    # OMIM and Orphanet IDs parsed from ClinVar's CLNDISDB field.
    diseases: frozenset[str]

    @staticmethod
    def of_json(record: _SpikeRecord) -> "SpikeVariant":
        return SpikeVariant(
            allele=Allele(record["chrom"], record["pos"], record["ref"], record["alt"]),
            clinvar_id=record["id"],
            gene=record["gene"],
            kind=SpikeKind(record["kind"]),
            diseases=frozenset(record["diseases"]),
        )

    def to_json(self) -> _SpikeRecord:
        return {
            "chrom": self.allele.chromosome,
            "pos": self.allele.position,
            "id": self.clinvar_id,
            "ref": self.allele.reference,
            "alt": self.allele.alternate,
            "gene": self.gene,
            "kind": self.kind.value,
            "diseases": sorted(self.diseases),
        }


@dataclass(frozen=True)
class Case:
    """A simulated patient's causal genotype and exact or noisy HPO terms."""

    name: str
    genotype: Genotype
    gene: str
    variants: tuple[SpikeVariant, ...]
    disease: str
    is_disease_named_by_clinvar: bool
    exact_terms: tuple[str, ...]
    noisy_terms: tuple[str, ...]

    @property
    def kinds(self) -> str:
        return "+".join(sorted(variant.kind.value for variant in self.variants))

    @staticmethod
    def of_json(record: _CaseRecord) -> "Case":
        return Case(
            name=record["case"],
            genotype=Genotype(record["type"]),
            gene=record["gene"],
            variants=tuple(map(SpikeVariant.of_json, record["variants"])),
            disease=record["disease"],
            is_disease_named_by_clinvar=record["disease_matches_clinvar"],
            exact_terms=tuple(record["exact"]),
            noisy_terms=tuple(record["noisy"]),
        )

    def to_json(self) -> _CaseRecord:
        return {
            "case": self.name,
            "type": self.genotype.value,
            "gene": self.gene,
            "variants": [variant.to_json() for variant in self.variants],
            "disease": self.disease,
            "disease_matches_clinvar": self.is_disease_named_by_clinvar,
            "exact": list(self.exact_terms),
            "noisy": list(self.noisy_terms),
        }


def read_cases(cases_path: Path) -> list[Case]:
    records: list[_CaseRecord] = json.loads(cases_path.read_text(encoding="utf8"))
    return list(map(Case.of_json, records))


def write_cases(cases: Sequence[Case], cases_path: Path) -> None:
    cases_path.write_text(
        json.dumps([case.to_json() for case in cases], indent=1), encoding="utf8"
    )


def spike_variant_of_clinvar_line(line: str) -> SpikeVariant | None:
    """Parse a pathogenic autosomal allele assigned to one gene, or return None.

    Require at least MINIMUM_REVIEW_STARS and an A/C/G/T allele within the length limit.
    """
    (
        chromosome,
        position,
        clinvar_id,
        reference,
        alternate,
        gene_field,
        consequence,
        significance,
        review_status,
        disease_field,
    ) = line.split("\t")
    genes = [entry.split(":")[0] for entry in gene_field.split("|") if entry != "."]
    is_usable = (
        chromosome in AUTOSOMES
        and significance in PATHOGENIC
        and parse_review_stars(review_status) >= MINIMUM_REVIEW_STARS
        and len(genes) == 1
        and set(reference + alternate) <= set("ACGT")
        and max(len(reference), len(alternate)) <= MAXIMUM_ALLELE_LENGTH
    )
    if not is_usable:
        return None
    if any(term in consequence for term in LOSS_OF_FUNCTION_CONSEQUENCES):
        kind = SpikeKind.LOSS_OF_FUNCTION
    elif "missense_variant" in consequence:
        kind = SpikeKind.MISSENSE
    else:
        kind = SpikeKind.OTHER
    diseases = set()
    for identifier in re.split(r"[,|]", disease_field):
        if identifier.startswith("OMIM:"):
            diseases.add(identifier)
        elif identifier.startswith("Orphanet:"):
            diseases.add("ORPHA:" + identifier.removeprefix("Orphanet:"))
    return SpikeVariant(
        allele=Allele(chromosome, int(position), reference, alternate),
        clinvar_id=clinvar_id,
        gene=genes[0],
        kind=kind,
        diseases=frozenset(diseases),
    )


def load_clinvar_pool(clinvar_path: Path) -> dict[str, list[SpikeVariant]]:
    lines = capture(
        [
            "bcftools", "query", "-i", 'INFO/CLNSIG~"athogenic"',
            "-f", CLINVAR_QUERY_FORMAT, str(clinvar_path),
        ]
    ).splitlines()  # fmt: skip
    pool: dict[str, list[SpikeVariant]] = defaultdict(list)
    for line in lines:
        variant = spike_variant_of_clinvar_line(line)
        if variant is not None:
            pool[variant.gene].append(variant)
    return pool


def is_phenotypic_abnormality(ontology: Ontology, term: str) -> bool:
    return (
        term != PHENOTYPIC_ABNORMALITY
        and PHENOTYPIC_ABNORMALITY in ontology.ancestors(term)
    )


def _sampleable_terms_of_disease(
    model: PhenotypeModel, gene: str
) -> dict[str, list[str]]:
    terms_of_disease = {}
    for disease, terms in model.disease_terms.get(gene, {}).items():
        abnormalities = sorted(
            filter(functools.partial(is_phenotypic_abnormality, model.ontology), terms)
        )
        if len(abnormalities) >= MINIMUM_TERMS_PER_DISEASE:
            terms_of_disease[disease] = abnormalities
    return terms_of_disease


def _make_case(
    genotype: Genotype,
    gene: str,
    variants: Sequence[SpikeVariant],
    *,
    model: PhenotypeModel,
    gene_disease_pairs: Sequence[tuple[str, str]],
    random_generator: random.Random,
) -> Case:
    ontology = model.ontology
    terms_of_disease = _sampleable_terms_of_disease(model, gene)
    clinvar_diseases = frozenset().union(*(variant.diseases for variant in variants))
    named_by_clinvar = sorted(clinvar_diseases & terms_of_disease.keys())
    disease = random_generator.choice(named_by_clinvar or sorted(terms_of_disease))
    exact_terms = random_generator.sample(terms_of_disease[disease], TERMS_PER_CASE)

    noisy_terms = []
    for term in exact_terms:
        parents = sorted(
            filter(
                functools.partial(is_phenotypic_abnormality, ontology),
                ontology.parents.get(term, ()),
            )
        )
        if parents and random_generator.random() < PARENT_TERM_PROBABILITY:
            noisy_terms.append(random_generator.choice(parents))
        else:
            noisy_terms.append(term)
    for _ in range(UNRELATED_TERMS_PER_CASE):
        other_gene, other_disease = random_generator.choice(gene_disease_pairs)
        other_terms = sorted(model.disease_terms[other_gene][other_disease])
        noisy_terms.append(random_generator.choice(other_terms))

    return Case(
        name=f"{genotype.value}_{gene}",
        genotype=genotype,
        gene=gene,
        variants=tuple(variants),
        disease=disease,
        is_disease_named_by_clinvar=bool(named_by_clinvar),
        exact_terms=tuple(exact_terms),
        noisy_terms=tuple(filter(ontology.is_feature, dict.fromkeys(noisy_terms))),
    )


def make_cases(
    pool: Mapping[str, Sequence[SpikeVariant]],
    gene_validity: Mapping[str, GeneValidity],
    model: PhenotypeModel,
    background_alleles: Collection[Allele],
    *,
    genes_per_group: int,
    random_generator: random.Random,
) -> list[Case]:
    """Sample recessive pairs, recessive homozygotes and dominant heterozygotes.

    Exclude background alleles and require genes with enough disease HPO terms. Each
    sampled recessive gene supplies a paired case and a homozygous case.
    """

    def usable_variants(gene: str) -> list[SpikeVariant]:
        return [
            variant for variant in pool[gene] if variant.allele not in background_alleles
        ]

    def validity(gene: str, inheritance: Inheritance) -> Validity:
        return gene_validity.get(gene, UNKNOWN_GENE).for_inheritance(inheritance)

    genes = [gene for gene in pool if _sampleable_terms_of_disease(model, gene)]
    recessive_genes = sorted(
        gene
        for gene in genes
        if validity(gene, Inheritance.AUTOSOMAL_RECESSIVE) >= Validity.STRONG
        and len({variant.allele.position for variant in usable_variants(gene)}) >= 2
    )
    dominant_genes = sorted(
        gene
        for gene in genes
        if validity(gene, Inheritance.AUTOSOMAL_DOMINANT) >= Validity.STRONG
        and validity(gene, Inheritance.AUTOSOMAL_RECESSIVE) < Validity.SUPPORTIVE
        and usable_variants(gene)
    )
    gene_disease_pairs = sorted(
        (gene, disease)
        for gene, terms_of_disease in model.disease_terms.items()
        for disease in terms_of_disease
    )

    make_case = functools.partial(
        _make_case,
        model=model,
        gene_disease_pairs=gene_disease_pairs,
        random_generator=random_generator,
    )
    cases = []
    for gene in random_generator.sample(recessive_genes, genes_per_group):
        variants = usable_variants(gene)
        first = random_generator.choice(variants)
        second = random_generator.choice(
            [
                variant
                for variant in variants
                if variant.allele.position != first.allele.position
            ]
        )
        cases.append(make_case(Genotype.COMPOUND_HETEROZYGOUS, gene, [first, second]))
        homozygous = random_generator.choice(variants)
        cases.append(make_case(Genotype.HOMOZYGOUS, gene, [homozygous]))
    for gene in random_generator.sample(dominant_genes, genes_per_group):
        heterozygous = random_generator.choice(usable_variants(gene))
        cases.append(make_case(Genotype.HETEROZYGOUS, gene, [heterozygous]))
    return cases
