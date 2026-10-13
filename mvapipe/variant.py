from dataclasses import dataclass
from enum import Enum
from typing import NamedTuple


class Allele(NamedTuple):
    chromosome: str
    position: int
    reference: str
    alternate: str

    def to_string(self) -> str:
        return f"{self.chromosome}:{self.position} {self.reference}>{self.alternate}"


class Zygosity(Enum):
    Heterozygous = "het"
    Homozygous = "hom"
    Hemizygous = "hemi"
    Other = "other"

    @staticmethod
    def of_genotype(genotype: str) -> "Zygosity":
        calls = alleles_of_genotype(genotype)
        alternates = [call for call in calls if call not in ("0", ".", "")]
        if len(calls) == 1:
            return Zygosity.Hemizygous if alternates else Zygosity.Other
        if len(alternates) == len(calls) and len(set(alternates)) == 1:
            return Zygosity.Homozygous
        if len(alternates) == 1:
            return Zygosity.Heterozygous
        return Zygosity.Other


def alleles_of_genotype(genotype: str) -> list[str]:
    return genotype.replace("|", "/").split("/")


@dataclass(frozen=True)
class Variant:
    """An allele in the sample, annotated for one gene."""

    allele: Allele
    gene: str
    consequence: str
    vep_impact: str
    exon: str
    genotype: str
    filter_status: str
    genotype_quality: float | None
    depth: float | None
    allele_depths: tuple[int, ...] | None
    phase_set: str | None
    phase_haplotype: int | None
    # CHROM|POS|REF|ALTs of the multi-allelic site the record was split from.
    original_site: str | None
    allele_frequency: float
    alphamissense: float | None
    clinvar_significance: str
    clinvar_review_status: str
    clinvar_conflict: str
    identifier: str
    is_gene_from_clinvar: bool
