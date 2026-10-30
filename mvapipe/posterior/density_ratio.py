import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from sklearn.linear_model import LogisticRegression

from mvapipe.posterior.features import find_clinvar_columns, remove_clinvar

MAXIMUM_ITERATIONS = 5000


@dataclass(frozen=True)
class DensityRatio:
    """A linear log density ratio: log(f1(x) / f0(x)) = x @ coefficients + intercept."""

    coefficients: NDArray[np.float64]
    intercept: float

    def log_ratio(self, features: NDArray[np.float64]) -> NDArray[np.float64]:
        return features @ self.coefficients + self.intercept


def fit_density_ratio(
    *,
    positive: NDArray[np.float64],
    negative: NDArray[np.float64],
    regularisation: float,
) -> DensityRatio:
    """Fit a log density ratio between positive and negative feature rows.

    Use `regularisation` as logistic regression's inverse penalty strength C. Correct
    the intercept by log(n_negative / n_positive) to remove the sample-count imbalance.
    """
    classifier = LogisticRegression(C=regularisation, max_iter=MAXIMUM_ITERATIONS)
    classifier.fit(
        np.vstack([positive, negative]),
        np.concatenate([np.ones(len(positive)), np.zeros(len(negative))]),
    )
    return DensityRatio(
        coefficients=classifier.coef_[0],
        intercept=float(classifier.intercept_[0])
        + math.log(len(negative) / len(positive)),
    )


@dataclass(frozen=True)
class StratifiedRatio:
    """Separate density ratios for alleles present in and absent from ClinVar.

    Add log((1 - p) / (1 - q)) for listed alleles and log(p / q) for unlisted alleles.
    Here p is the prior probability that a causal allele is unlisted, and q is the
    fraction of background alleles that are unlisted.
    """

    listed: DensityRatio
    unlisted: DensityRatio
    log_listed_share: float
    log_unlisted_share: float
    absent_column: int

    def log_ratio(self, features: NDArray[np.float64]) -> NDArray[np.float64]:
        is_unlisted = features[:, self.absent_column] == 1
        return np.where(
            is_unlisted,
            self.unlisted.log_ratio(features) + self.log_unlisted_share,
            self.listed.log_ratio(features) + self.log_listed_share,
        )


def fit_stratified_ratio(
    *,
    positive: NDArray[np.float64],
    negative: NDArray[np.float64],
    frequency_bounds: tuple[float, ...],
    listed_regularisation: float,
    unlisted_regularisation: float,
    unlisted_prior: float,
) -> StratifiedRatio:
    """Fit separate ratios for alleles present in and absent from ClinVar.

    Compare listed causal alleles with listed background alleles. For the unlisted fit,
    remove ClinVar features from all causal alleles and compare them with unlisted
    background alleles.
    """
    _, absent_column = find_clinvar_columns(frequency_bounds)
    is_listed = positive[:, absent_column] == 0
    is_null_listed = negative[:, absent_column] == 0
    unlisted_share = float(np.mean(~is_null_listed))
    return StratifiedRatio(
        listed=fit_density_ratio(
            positive=positive[is_listed],
            negative=negative[is_null_listed],
            regularisation=listed_regularisation,
        ),
        unlisted=fit_density_ratio(
            positive=remove_clinvar(positive, frequency_bounds),
            negative=negative[~is_null_listed],
            regularisation=unlisted_regularisation,
        ),
        log_listed_share=math.log((1 - unlisted_prior) / (1 - unlisted_share)),
        log_unlisted_share=math.log(unlisted_prior / unlisted_share),
        absent_column=absent_column,
    )
