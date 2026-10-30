import math
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum

from mvapipe.gene_validity import UNKNOWN_GENE, GeneValidity, Inheritance, Validity

PSEUDO_COUNT = 0.5


class Genotype(Enum):
    DOMINANT = "dominant"
    HOMOZYGOUS = "homozygous"
    COMPOUND = "compound"
    # One small allele, with a second causal allele absent from the small-variant records.
    SINGLE_RECESSIVE = "single_recessive"


RECESSIVE_GENOTYPES = (Genotype.HOMOZYGOUS, Genotype.COMPOUND, Genotype.SINGLE_RECESSIVE)


class ValidityClass(Enum):
    STRONG = "strong"
    MODERATE = "moderate"
    LIMITED = "limited"
    NONE = "none"


def classify_validity(validity: GeneValidity, inheritance: Inheritance) -> ValidityClass:
    """Group a gene's GenCC classification for the requested inheritance mode."""
    level = validity.for_inheritance(inheritance)
    if level >= Validity.STRONG:
        return ValidityClass.STRONG
    if level == Validity.MODERATE:
        return ValidityClass.MODERATE
    if level >= Validity.LIMITED:
        return ValidityClass.LIMITED
    return ValidityClass.NONE


@dataclass(frozen=True)
class Priors:
    log_mode: dict[Inheritance, float]
    log_genotype: dict[Genotype, float]
    # Per-gene log P(gene | mode), shared by GenCC class and normalised over all genes.
    log_gene_by_class: dict[Inheritance, dict[ValidityClass, float]]
    log_unexplained: float
    log_explained: float

    def log_gene_mode(
        self,
        gene: str,
        inheritance: Inheritance,
        gene_validity: Mapping[str, GeneValidity],
    ) -> float:
        """Return log P(explained, mode, gene)."""
        validity = gene_validity.get(gene, UNKNOWN_GENE)
        validity_class = classify_validity(validity, inheritance)
        return (
            self.log_explained
            + self.log_mode[inheritance]
            + self.log_gene_by_class[inheritance][validity_class]
        )


def estimate_priors(
    training: Sequence[tuple[str, Genotype]],
    gene_universe: Sequence[str],
    gene_validity: Mapping[str, GeneValidity],
    *,
    unexplained_prior: float,
) -> Priors:
    """Estimate inheritance, genotype and gene priors from the training cases.

    `unexplained_prior` reserves probability for cases unexplained by the genome's
    variants. Normalise gene priors over `gene_universe`.
    """
    modes = Counter(
        Inheritance.AUTOSOMAL_RECESSIVE
        if genotype in RECESSIVE_GENOTYPES
        else Inheritance.AUTOSOMAL_DOMINANT
        for _, genotype in training
    )
    genotypes = Counter(
        genotype for _, genotype in training if genotype in RECESSIVE_GENOTYPES
    )
    recessive_count = sum(genotypes.values())
    log_gene_by_class = {}
    for inheritance in (Inheritance.AUTOSOMAL_DOMINANT, Inheritance.AUTOSOMAL_RECESSIVE):
        is_recessive = inheritance is Inheritance.AUTOSOMAL_RECESSIVE
        genes_in_class = Counter(
            classify_validity(gene_validity.get(gene, UNKNOWN_GENE), inheritance)
            for gene in gene_universe
        )
        cases_in_class = Counter(
            classify_validity(gene_validity.get(gene, UNKNOWN_GENE), inheritance)
            for gene, genotype in training
            if (genotype in RECESSIVE_GENOTYPES) == is_recessive
        )
        weights = {
            name: (cases_in_class[name] + PSEUDO_COUNT) / max(genes_in_class[name], 1)
            for name in ValidityClass
        }
        normaliser = sum(weights[name] * genes_in_class[name] for name in ValidityClass)
        log_gene_by_class[inheritance] = {
            name: math.log(weights[name] / normaliser) for name in ValidityClass
        }
    case_count = sum(modes.values())
    return Priors(
        log_mode={mode: math.log(count / case_count) for mode, count in modes.items()},
        log_genotype={
            genotype: math.log(
                (genotypes[genotype] + PSEUDO_COUNT)
                / (recessive_count + PSEUDO_COUNT * len(RECESSIVE_GENOTYPES))
            )
            for genotype in RECESSIVE_GENOTYPES
        },
        log_gene_by_class=log_gene_by_class,
        log_unexplained=math.log(unexplained_prior),
        log_explained=math.log(1 - unexplained_prior),
    )
