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
    HETEROZYGOUS = "het"
    HOMOZYGOUS = "hom"
    HEMIZYGOUS = "hemi"
    OTHER = "other"

    @staticmethod
    def from_genotype(genotype: str) -> "Zygosity":
        """Classify a VCF genotype such as '0/1', '1|1' or '1'."""
        calls = split_genotype(genotype)
        if any(call in (".", "") for call in calls) or len(calls) > 2:
            return Zygosity.OTHER
        alternates = [call for call in calls if call not in ("0", ".", "")]
        if len(calls) == 1:
            return Zygosity.HEMIZYGOUS if alternates else Zygosity.OTHER
        if len(alternates) == len(calls) and len(set(alternates)) == 1:
            return Zygosity.HOMOZYGOUS
        if len(alternates) == 1:
            return Zygosity.HETEROZYGOUS
        return Zygosity.OTHER


def split_genotype(genotype: str) -> list[str]:
    """Return the allele indices of a phased or unphased VCF genotype."""
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
    # Original site before splitting: CHROM|POS|REF|ALTs, without the allele index.
    original_site: str | None
    allele_frequency: float | None
    alphamissense: float | None
    clinvar_significance: str
    clinvar_review_status: str
    clinvar_conflict: str
    identifier: str
    is_gene_from_clinvar: bool
