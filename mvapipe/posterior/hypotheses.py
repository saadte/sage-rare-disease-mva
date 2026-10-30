import math
from collections import defaultdict
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.special import logsumexp

from mvapipe.gene_validity import GeneValidity, Inheritance
from mvapipe.posterior.density_ratio import StratifiedRatio
from mvapipe.posterior.phenotype import GenePhenotype
from mvapipe.posterior.priors import Genotype, Priors
from mvapipe.posterior.training import Context, GeneShares
from mvapipe.ranking import MITOCHONDRIAL, Phase, infer_phase
from mvapipe.variant import Allele, Variant, Zygosity

MAXIMUM_HETEROZYGOTES_PAIRED = 30
MAXIMUM_SHAPE = 1000.0
# An empty allele set represents a case unexplained by variants in this genome.
UNEXPLAINED: frozenset[Allele] = frozenset()


@dataclass(frozen=True)
class Genome:
    """Nuclear variants with zygosity masks that group hemizygotes with homozygotes."""

    variants: tuple[Variant, ...]
    is_heterozygous: NDArray[np.bool_]
    is_homozygous: NDArray[np.bool_]
    log_quality: NDArray[np.float64] | None = None

    def __post_init__(self) -> None:
        count = len(self.variants)
        if (
            self.is_heterozygous.shape != (count,)
            or self.is_homozygous.shape != (count,)
            or (self.log_quality is not None and self.log_quality.shape != (count,))
        ):
            message = "genotype masks and quality values must align with genome variants"
            raise ValueError(message)

    @staticmethod
    def from_variants(variants: Sequence[Variant]) -> "Genome":
        kept = tuple(
            variant
            for variant in variants
            if variant.allele.chromosome not in MITOCHONDRIAL
        )
        zygosity = [Zygosity.from_genotype(variant.genotype) for variant in kept]
        return Genome(
            variants=kept,
            is_heterozygous=np.array(
                [z is Zygosity.HETEROZYGOUS for z in zygosity], dtype=bool
            ),
            is_homozygous=np.array(
                [z in (Zygosity.HOMOZYGOUS, Zygosity.HEMIZYGOUS) for z in zygosity],
                dtype=bool,
            ),
        )


@dataclass(frozen=True)
class BackgroundRates:
    """Gamma priors for each gene's heterozygous and homozygous background row rates.

    Set prior means from gnomAD gene shares and genome-wide row counts in those genes.
    Update each rate with the gene's observed rows, removing hypothesis alleles in turn.
    """

    gene_shares: GeneShares
    heterozygous_count: int
    homozygous_count: int
    heterozygous_shape: float
    homozygous_shape: float

    def expect(self, gene: str, *, is_homozygous: bool) -> float:
        """Return the gene's expected row count under the prior."""
        count = self.homozygous_count if is_homozygous else self.heterozygous_count
        return count * self.gene_shares.shares.get(gene, self.gene_shares.default_share)

    def compute_log_rate(
        self, gene: str, *, observed: int, allele_count: int, is_homozygous: bool
    ) -> float:
        """Sum log posterior mean rates while removing hypothesis alleles in turn.

        `observed` includes the hypothesis alleles. For allele j, update the gamma prior
        with observed - j rows, where j runs from 1 through `allele_count`.
        """
        mean = self.expect(gene, is_homozygous=is_homozygous)
        shape = self.homozygous_shape if is_homozygous else self.heterozygous_shape
        return sum(
            math.log((shape + observed - j) / (shape / mean + 1))
            for j in range(1, allele_count + 1)
        )


def _count_gene_rows(genome: Genome, *, is_homozygous: bool) -> dict[str, int]:
    is_row = genome.is_homozygous if is_homozygous else genome.is_heterozygous
    counts: dict[str, int] = defaultdict(int)
    for variant, is_counted in zip(genome.variants, is_row, strict=True):
        if is_counted:
            counts[variant.gene] += 1
    return dict(counts)


def _estimate_shape(
    counts: Mapping[str, int], gene_shares: GeneShares, row_count: int
) -> float:
    """Estimate the gamma shape by method of moments, capped at MAXIMUM_SHAPE."""
    excess, scale = 0.0, 0.0
    for gene, share in gene_shares.shares.items():
        mean = row_count * share
        observed = counts.get(gene, 0)
        excess += (observed - mean) ** 2 - observed
        scale += mean**2
    return 1 / max(excess / scale, 1 / MAXIMUM_SHAPE)


def estimate_background_rates(genome: Genome, gene_shares: GeneShares) -> BackgroundRates:
    """Estimate background row rates using genes with gnomAD shares.

    Raise ValueError if these genes lack heterozygous or homozygous rows.
    """
    is_in_universe = np.array(
        [variant.gene in gene_shares.shares for variant in genome.variants], dtype=bool
    )
    heterozygous_count = int((genome.is_heterozygous & is_in_universe).sum())
    homozygous_count = int((genome.is_homozygous & is_in_universe).sum())
    if not heterozygous_count or not homozygous_count:
        message = (
            "the posterior method needs heterozygous and homozygous variants "
            "in genes with a gnomAD share"
        )
        raise ValueError(message)
    return BackgroundRates(
        gene_shares=gene_shares,
        heterozygous_count=heterozygous_count,
        homozygous_count=homozygous_count,
        heterozygous_shape=_estimate_shape(
            _count_gene_rows(genome, is_homozygous=False), gene_shares, heterozygous_count
        ),
        homozygous_shape=_estimate_shape(
            _count_gene_rows(genome, is_homozygous=True), gene_shares, homozygous_count
        ),
    )


@dataclass(frozen=True)
class Hypotheses:
    """Allele sets and genes with log scores before phenotype evidence is added."""

    log_score: NDArray[np.float64]
    genes: list[str]
    allele_sets: list[frozenset[Allele]]


def _pair_heterozygotes(
    genome: Genome, heterozygous: Sequence[int], log_recessive: NDArray[np.float64]
) -> Iterator[tuple[int, int, float]]:
    """Yield distinct heterozygous row pairs and their log phase factors.

    Pair at most MAXIMUM_HETEROZYGOTES_PAIRED rows with the highest recessive log ratios.
    Skip cis pairs; use log(2) for trans pairs and 0 for unknown phase.
    """
    paired = sorted(heterozygous, key=lambda row: -log_recessive[row])
    paired = paired[:MAXIMUM_HETEROZYGOTES_PAIRED]
    for index, first in enumerate(paired):
        for second in paired[index + 1 :]:
            one, two = genome.variants[first], genome.variants[second]
            phase = infer_phase(one, two)
            if one.allele == two.allele or phase is Phase.CIS:
                continue
            yield first, second, math.log(2) if phase is Phase.TRANS else 0.0


def _compute_gene_log_rates(
    rates: BackgroundRates, gene: str, *, heterozygous_count: int, homozygous_count: int
) -> tuple[float, float, float]:
    """Return log background rates for a heterozygote, a pair and a homozygote.

    Use a placeholder log rate of 0 when too few rows exist for that hypothesis.
    """
    log_rate_heterozygous = log_rate_pair = log_rate_homozygous = 0.0
    if heterozygous_count:
        log_rate_heterozygous = rates.compute_log_rate(
            gene, observed=heterozygous_count, allele_count=1, is_homozygous=False
        )
    if heterozygous_count >= 2:
        log_rate_pair = rates.compute_log_rate(
            gene, observed=heterozygous_count, allele_count=2, is_homozygous=False
        )
    if homozygous_count:
        log_rate_homozygous = rates.compute_log_rate(
            gene, observed=homozygous_count, allele_count=1, is_homozygous=True
        )
    return log_rate_heterozygous, log_rate_pair, log_rate_homozygous


def enumerate_hypotheses(
    genome: Genome,
    features: NDArray[np.float64],
    ratios: Mapping[Context, StratifiedRatio],
    rates: BackgroundRates,
    priors: Priors,
    gene_validity: Mapping[str, GeneValidity],
) -> Hypotheses:
    """Build dominant, single recessive, homozygous and compound hypotheses per gene.

    Use `_pair_heterozygotes` for compound pairs. Expect feature rows in genome variant
    order; raise ValueError if the row count differs from the variant count.
    """
    if len(features) != len(genome.variants):
        message = f"{len(features)} feature rows for {len(genome.variants)} variants"
        raise ValueError(message)
    quality = 0.0 if genome.log_quality is None else genome.log_quality
    log_dominant = ratios[Context.DOMINANT].log_ratio(features) + quality
    log_recessive = ratios[Context.RECESSIVE_HETEROZYGOUS].log_ratio(features) + quality
    log_homozygous = ratios[Context.RECESSIVE_HOMOZYGOUS].log_ratio(features) + quality
    rows_of_gene: dict[str, list[int]] = defaultdict(list)
    for row, variant in enumerate(genome.variants):
        rows_of_gene[variant.gene].append(row)
    scores: list[float] = []
    genes: list[str] = []
    allele_sets: list[frozenset[Allele]] = []

    def add(gene: str, score: float, alleles: frozenset[Allele]) -> None:
        scores.append(score)
        genes.append(gene)
        allele_sets.append(alleles)

    for gene, rows in rows_of_gene.items():
        dominant = priors.log_gene_mode(
            gene, Inheritance.AUTOSOMAL_DOMINANT, gene_validity
        )
        recessive = priors.log_gene_mode(
            gene, Inheritance.AUTOSOMAL_RECESSIVE, gene_validity
        )
        heterozygous = [row for row in rows if genome.is_heterozygous[row]]
        homozygous_rows = [row for row in rows if genome.is_homozygous[row]]
        log_rate_heterozygous, log_rate_pair, log_rate_homozygous = (
            _compute_gene_log_rates(
                rates,
                gene,
                heterozygous_count=len(heterozygous),
                homozygous_count=len(homozygous_rows),
            )
        )
        single = recessive + priors.log_genotype[Genotype.SINGLE_RECESSIVE]
        homozygous = recessive + priors.log_genotype[Genotype.HOMOZYGOUS]
        compound = recessive + priors.log_genotype[Genotype.COMPOUND]
        for row in heterozygous:
            alleles = frozenset({genome.variants[row].allele})
            add(gene, dominant + log_dominant[row] - log_rate_heterozygous, alleles)
            add(gene, single + log_recessive[row] - log_rate_heterozygous, alleles)
        for row in homozygous_rows:
            alleles = frozenset({genome.variants[row].allele})
            add(gene, homozygous + log_homozygous[row] - log_rate_homozygous, alleles)
        pairs = _pair_heterozygotes(genome, heterozygous, log_recessive)
        for first, second, log_phase in pairs:
            score = (
                compound + log_recessive[first] + log_recessive[second]
                - log_rate_pair + log_phase
            )  # fmt: skip
            alleles = frozenset(
                {genome.variants[first].allele, genome.variants[second].allele}
            )
            add(gene, score, alleles)
    return Hypotheses(log_score=np.array(scores), genes=genes, allele_sets=allele_sets)


@dataclass(frozen=True)
class RankedAlleleSets:
    """Allele sets in descending posterior order, plus the unexplained probability."""

    allele_sets: list[frozenset[Allele]]
    posteriors: NDArray[np.float64]
    unexplained: float
    gene_posteriors: dict[str, float] = field(default_factory=dict)
    is_complete: bool = True


def rank_allele_sets(
    hypotheses: Hypotheses,
    gene_phenotype: GenePhenotype | None,
    priors: Priors,
    *,
    top: int | None,
) -> RankedAlleleSets:
    """Sum hypothesis posteriors by allele set and return the `top` sets.

    Include the unexplained hypothesis when normalising. If `gene_phenotype` is None,
    use a phenotype likelihood ratio of 1 for every gene.
    """
    log_score = hypotheses.log_score
    if gene_phenotype is not None:
        log_score = log_score + np.array(
            [gene_phenotype.log_ratio(gene) for gene in hypotheses.genes]
        )
    normaliser = logsumexp(np.append(log_score, priors.log_unexplained))
    probability = np.exp(log_score - normaliser)
    totals: dict[frozenset[Allele], float] = defaultdict(float)
    gene_totals: dict[str, float] = defaultdict(float)
    for index in np.argsort(-probability, kind="stable"):
        posterior = float(probability[index])
        totals[hypotheses.allele_sets[index]] += posterior
        gene_totals[hypotheses.genes[index]] += posterior
    ordered = sorted(totals.items(), key=lambda item: -item[1])[:top]
    return RankedAlleleSets(
        allele_sets=[allele_set for allele_set, _ in ordered],
        posteriors=np.array([posterior for _, posterior in ordered]),
        unexplained=float(np.exp(priors.log_unexplained - normaliser)),
        gene_posteriors=dict(gene_totals),
        is_complete=top is None or top >= len(totals),
    )
