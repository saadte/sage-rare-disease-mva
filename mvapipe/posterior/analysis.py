import csv
import random
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from mvapipe.posterior.density_ratio import StratifiedRatio, fit_stratified_ratio
from mvapipe.posterior.features import FREQUENCY_SCHEMES, build_feature_matrix
from mvapipe.posterior.folds import assign_gene_folds
from mvapipe.posterior.hypotheses import (
    UNEXPLAINED,
    Genome,
    RankedAlleleSets,
    enumerate_hypotheses,
    estimate_background_rates,
    rank_allele_sets,
)
from mvapipe.posterior.phenotype import fit_phenotype_ratio, score_patient_phenotype
from mvapipe.posterior.priors import estimate_priors
from mvapipe.posterior.quality import compute_log_quality_ratio
from mvapipe.posterior.selection import (
    PhenotypeSettings,
    Stratum,
    VariantSettings,
    select_phenotype_settings,
    select_variant_settings,
)
from mvapipe.posterior.training import (
    SEED,
    Context,
    TrainingData,
    build_phenotype_table,
    collect_positive_variants,
)
from mvapipe.submission import Prediction
from mvapipe.variant import Allele, Variant

UNEXPLAINED_PRIOR = 0.5
# Prior probability that a causal allele is absent from ClinVar.
UNLISTED_PRIOR = 0.75


@dataclass(frozen=True)
class PosteriorConfig:
    """Base weights and optional quality evidence for the fitted ranking method."""

    unexplained_prior: float = UNEXPLAINED_PRIOR
    unlisted_prior: float = UNLISTED_PRIOR
    use_quality: bool = False

    def __post_init__(self) -> None:
        for name, value in (
            ("unexplained_prior", self.unexplained_prior),
            ("unlisted_prior", self.unlisted_prior),
        ):
            if not 0 < value < 1:
                message = f"{name} must be strictly between 0 and 1"
                raise ValueError(message)


DEFAULT_CONFIG = PosteriorConfig()


@dataclass(frozen=True)
class PosteriorAnalysis:
    variant_settings: VariantSettings
    # None when phenotype scoring is disabled.
    phenotype_settings: PhenotypeSettings | None
    ranked: RankedAlleleSets
    genes_of_allele: dict[Allele, frozenset[str]]


def analyse_genome(
    variants: Sequence[Variant],
    patient_terms: Sequence[str] | None,
    training: TrainingData,
    *,
    top: int | None,
    config: PosteriorConfig = DEFAULT_CONFIG,
    patient_publications: frozenset[str] | None = None,
) -> PosteriorAnalysis:
    """Fit posterior models using the patient's genome as background, then rank alleles.

    Choose settings by cross-validation with genes kept in separate folds. Omit the
    phenotype likelihood ratio when `patient_terms` is None. For literature evaluation,
    pass `patient_publications` to exclude annotations sourced only from that case.
    Raise ValueError if
    background rates cannot be estimated or a setting has fewer than two scored folds.
    """
    genome = Genome.from_variants(variants)
    if config.use_quality:
        genome = replace(genome, log_quality=compute_log_quality_ratio(genome))
    rates = estimate_background_rates(genome, training.gene_shares)
    gene_folds = assign_gene_folds(
        [case.gene for case in training.cases]
        + [variant.gene for variant in genome.variants],
        known_folds=training.gene_folds,
    )
    folds = sorted({gene_folds[case.gene] for case in training.cases})
    positives = collect_positive_variants(training.cases, training.alleles, gene_folds)
    variant_settings = select_variant_settings(
        positives,
        genome.variants,
        null_folds=np.array([gene_folds[variant.gene] for variant in genome.variants]),
        null_is_heterozygous=genome.is_heterozygous,
        null_is_homozygous=genome.is_homozygous,
        folds=folds,
    )

    phenotype_settings = None
    gene_phenotype = None
    if patient_terms is not None:
        table = build_phenotype_table(training, gene_folds, random.Random(SEED))
        phenotype_settings = select_phenotype_settings(table, folds)
        phenotype_ratio = fit_phenotype_ratio(
            table,
            similarity_bins=phenotype_settings.similarity_bins,
            regularisation=phenotype_settings.regularisation,
        )
        gene_phenotype = score_patient_phenotype(
            phenotype_ratio,
            training.similarity,
            patient_terms,
            leave_source_out=training.leave_source_out
            if patient_publications is not None
            else None,
            publications=patient_publications or frozenset(),
        )

    bounds = FREQUENCY_SCHEMES[variant_settings.frequency_scheme]
    features = build_feature_matrix(genome.variants, bounds)
    ratios: dict[Context, StratifiedRatio] = {}
    for context in Context:
        is_null_row = (
            genome.is_homozygous
            if context is Context.RECESSIVE_HOMOZYGOUS
            else genome.is_heterozygous
        )
        ratios[context] = fit_stratified_ratio(
            positive=build_feature_matrix(
                [variant for variant, _ in positives[context]], bounds
            ),
            negative=features[is_null_row],
            frequency_bounds=bounds,
            listed_regularisation=variant_settings.regularisation[
                context, Stratum.LISTED
            ],
            unlisted_regularisation=variant_settings.regularisation[
                context, Stratum.UNLISTED
            ],
            unlisted_prior=config.unlisted_prior,
        )
    priors = estimate_priors(
        [(case.gene, case.genotype) for case in training.cases],
        training.gene_universe,
        training.gene_validity,
        unexplained_prior=config.unexplained_prior,
    )
    hypotheses = enumerate_hypotheses(
        genome,
        features,
        ratios,
        rates,
        priors,
        training.gene_validity,
    )
    genes_of_allele: dict[Allele, set[str]] = defaultdict(set)
    for variant in genome.variants:
        genes_of_allele[variant.allele].add(variant.gene)
    return PosteriorAnalysis(
        variant_settings=variant_settings,
        phenotype_settings=phenotype_settings,
        ranked=rank_allele_sets(hypotheses, gene_phenotype, priors, top=top),
        genes_of_allele={
            allele: frozenset(genes) for allele, genes in genes_of_allele.items()
        },
    )


def _find_genes(analysis: PosteriorAnalysis, alleles: frozenset[Allele]) -> list[str]:
    return sorted(
        frozenset().union(*(analysis.genes_of_allele[allele] for allele in alleles))
    )


def convert_allele_sets(analysis: PosteriorAnalysis) -> list[Prediction]:
    """Convert ranked allele sets to predictions, using posterior probability as EPCR."""
    ranked = analysis.ranked
    return [
        Prediction(
            alleles=alleles,
            score=posterior,
            epcr=posterior,
            note=" ".join(_find_genes(analysis, alleles)),
        )
        for alleles, posterior in zip(
            ranked.allele_sets, ranked.posteriors.tolist(), strict=True
        )
    ]


def write_posterior_table(
    analysis: PosteriorAnalysis, table_path: Path, *, limit: int | None = None
) -> None:
    """Write a TSV of allele sets and the unexplained hypothesis, highest posterior first.

    Leave the variant columns empty for the unexplained hypothesis.
    """
    ranked = analysis.ranked
    entries = [
        *zip(ranked.allele_sets[:limit], ranked.posteriors[:limit].tolist(), strict=True),
        (UNEXPLAINED, ranked.unexplained),
    ]
    entries.sort(key=lambda entry: -entry[1])
    with table_path.open("w", encoding="utf-8", newline="") as table:
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(["rank", "posterior", "genes", "variant_1", "variant_2"])
        for rank, (alleles, posterior) in enumerate(entries, start=1):
            genes = _find_genes(analysis, alleles)
            variants = [allele.to_string() for allele in sorted(alleles)]
            variants += [""] * (2 - len(variants))
            writer.writerow([rank, format(posterior, ".17g"), ",".join(genes), *variants])
