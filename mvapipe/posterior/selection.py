import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from typing import TypeVar

import numpy as np
from numpy.typing import NDArray
from scipy.special import logsumexp

from mvapipe.posterior.density_ratio import fit_density_ratio
from mvapipe.posterior.features import (
    FREQUENCY_SCHEMES,
    build_feature_matrix,
    find_clinvar_columns,
    remove_clinvar,
)
from mvapipe.posterior.phenotype import (
    SIMILARITY_SCHEMES,
    PhenotypeTable,
    build_similarity_features,
)
from mvapipe.posterior.training import Context
from mvapipe.variant import Variant

# Logistic regression C values; smaller values impose stronger regularisation.
VARIANT_GRID = (1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0)
PHENOTYPE_GRID = (0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0)
SELECTION_NEGATIVES_PER_CASE = 50

Key = TypeVar("Key")


class Stratum(Enum):
    """ClinVar listing status used to split variant model fits."""

    LISTED = "listed"
    UNLISTED = "unlisted"


@dataclass(frozen=True)
class VariantSettings:
    frequency_scheme: str
    regularisation: dict[tuple[Context, Stratum], float]


@dataclass(frozen=True)
class PhenotypeSettings:
    similarity_bins: int
    regularisation: float


def compute_heldout_log_likelihood(
    *, log_ratio_positive: NDArray[np.float64], log_ratio_negative: NDArray[np.float64]
) -> float:
    """Compute mean(log positive ratio) - log(mean negative ratio).

    Adding the same constant to all input log ratios leaves the score unchanged.
    """
    mean_log_ratio_negative = logsumexp(log_ratio_negative) - math.log(
        len(log_ratio_negative)
    )
    return float(log_ratio_positive.mean() - mean_log_ratio_negative)


def compute_cross_validated_scores(
    *,
    positive: NDArray[np.float64],
    positive_fold: NDArray[np.int_],
    negative: NDArray[np.float64],
    negative_fold: NDArray[np.int_],
    folds: Sequence[int],
    regularisation: float,
) -> NDArray[np.float64]:
    """Score held-out folds; use NaN if test rows or training positives are missing."""
    scores = np.full(len(folds), np.nan)
    for index, fold in enumerate(folds):
        test_positive = positive[positive_fold == fold]
        test_negative = negative[negative_fold == fold]
        train_positive = positive[positive_fold != fold]
        if not len(test_positive) or not len(test_negative) or not len(train_positive):
            continue
        ratio = fit_density_ratio(
            positive=train_positive,
            negative=negative[negative_fold != fold],
            regularisation=regularisation,
        )
        scores[index] = compute_heldout_log_likelihood(
            log_ratio_positive=ratio.log_ratio(test_positive),
            log_ratio_negative=ratio.log_ratio(test_negative),
        )
    return scores


def _is_within_one_standard_error(shortfall: NDArray[np.float64]) -> bool:
    present = shortfall[~np.isnan(shortfall)]
    standard_error = present.std(ddof=1) / math.sqrt(len(present))
    return bool(present.mean() <= standard_error)


def find_simplest_within_one_standard_error(
    scores: Mapping[Key, NDArray[np.float64]], order: Sequence[Key]
) -> Key:
    """Select the simplest setting within one standard error of the best.

    Supply keys in `order` from simplest to most complex. Compare each setting with the
    best mean score using paired fold differences: their mean shortfall must be no
    greater than their standard error. Raise ValueError if any setting has fewer than
    two scored folds.
    """
    if min(np.count_nonzero(~np.isnan(values)) for values in scores.values()) < 2:
        message = "cross-validation scored fewer than two folds of a setting"
        raise ValueError(message)
    best = max(scores, key=lambda key: np.nanmean(scores[key]))
    return next(
        key for key in order if _is_within_one_standard_error(scores[best] - scores[key])
    )


def select_variant_settings(
    positives: Mapping[Context, Sequence[tuple[Variant, int]]],
    null_variants: Sequence[Variant],
    *,
    null_folds: NDArray[np.int_],
    null_is_heterozygous: NDArray[np.bool_],
    null_is_homozygous: NDArray[np.bool_],
    folds: Sequence[int],
) -> VariantSettings:
    """Choose frequency bins and regularisation by gene-fold cross-validation.

    Fit each genotype context separately for listed and unlisted ClinVar alleles, as in
    `fit_stratified_ratio`. First choose the coarsest frequency scheme within one
    standard error of the best combined score, with each model at its best C. Then choose
    the strongest regularisation within one standard error for each model.
    """
    scores: dict[tuple[str, Context, Stratum], dict[float, NDArray[np.float64]]] = {}
    for scheme, bounds in FREQUENCY_SCHEMES.items():
        null_features = build_feature_matrix(null_variants, bounds)
        _, absent_column = find_clinvar_columns(bounds)
        is_null_listed = null_features[:, absent_column] == 0
        for context in Context:
            positive = build_feature_matrix(
                [variant for variant, _ in positives[context]], bounds
            )
            positive_fold = np.array([fold for _, fold in positives[context]])
            is_null_row = (
                null_is_homozygous
                if context is Context.RECESSIVE_HOMOZYGOUS
                else null_is_heterozygous
            )
            is_listed = positive[:, absent_column] == 0
            strata = [
                (
                    Stratum.LISTED,
                    positive[is_listed],
                    positive_fold[is_listed],
                    is_null_row & is_null_listed,
                ),
                (
                    Stratum.UNLISTED,
                    remove_clinvar(positive, bounds),
                    positive_fold,
                    is_null_row & ~is_null_listed,
                ),
            ]
            for stratum, model_positive, model_fold, is_model_row in strata:
                scores[scheme, context, stratum] = {
                    regularisation: compute_cross_validated_scores(
                        positive=model_positive,
                        positive_fold=model_fold,
                        negative=null_features[is_model_row],
                        negative_fold=null_folds[is_model_row],
                        folds=folds,
                        regularisation=regularisation,
                    )
                    for regularisation in VARIANT_GRID
                }
    models = [(context, stratum) for context in Context for stratum in Stratum]
    totals = {
        scheme: np.sum(
            [max(scores[scheme, *model].values(), key=np.nanmean) for model in models],
            axis=0,
        )
        for scheme in FREQUENCY_SCHEMES
    }
    scheme = find_simplest_within_one_standard_error(totals, list(FREQUENCY_SCHEMES))
    regularisation = {
        model: find_simplest_within_one_standard_error(
            scores[scheme, *model], VARIANT_GRID
        )
        for model in models
    }
    return VariantSettings(frequency_scheme=scheme, regularisation=regularisation)


def select_phenotype_settings(
    table: PhenotypeTable, folds: Sequence[int]
) -> PhenotypeSettings:
    """Choose similarity bins and regularisation within one standard error of the best.

    Prefer fewer bins, then stronger regularisation. Use each case's causal gene and
    its first SELECTION_NEGATIVES_PER_CASE sampled genes.
    """
    keep = table.is_causal | (table.order <= SELECTION_NEGATIVES_PER_CASE)
    is_causal = table.is_causal[keep]
    fold = table.fold[keep]
    scores: dict[tuple[int, float], NDArray[np.float64]] = {}
    for bins, bounds in SIMILARITY_SCHEMES.items():
        features = build_similarity_features(
            table.similarity[keep], table.is_annotated[keep], bounds
        )
        for regularisation in PHENOTYPE_GRID:
            scores[bins, regularisation] = compute_cross_validated_scores(
                positive=features[is_causal],
                positive_fold=fold[is_causal],
                negative=features[~is_causal],
                negative_fold=fold[~is_causal],
                folds=folds,
                regularisation=regularisation,
            )
    bins, regularisation = find_simplest_within_one_standard_error(scores, sorted(scores))
    return PhenotypeSettings(similarity_bins=bins, regularisation=regularisation)
