"""Optional call-quality evidence fitted against common variants in the same genome."""

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from mvapipe.posterior.density_ratio import fit_density_ratio
from mvapipe.posterior.hypotheses import Genome
from mvapipe.variant import Variant

BALANCE_BOUNDS = (0.2, 0.3, 0.7, 0.8)
DEPTH_BOUNDS = (8, 15, 30, 60)
QUALITY_BOUNDS = (20, 40, 60)
COMMON_FREQUENCY = 0.01
COLUMN_COUNT = (
    sum(len(bounds) + 2 for bounds in (BALANCE_BOUNDS, DEPTH_BOUNDS, QUALITY_BOUNDS)) + 1
)


def _one_hot(value: float | None, bounds: Sequence[float]) -> list[float]:
    row = [0.0] * (len(bounds) + 2)
    index = (
        len(row) - 1
        if value is None
        else int(np.searchsorted(bounds, value, side="right"))
    )
    row[index] = 1.0
    return row


def build_quality_features(variants: Sequence[Variant]) -> NDArray[np.float64]:
    """Encode allele balance, depth, genotype quality and filter status by variant.

    Each numeric variable has its own missing-value column. Use summed allele depths
    when DP is missing. Columns follow balance, depth, quality, then failed filter.
    """
    rows = []
    for variant in variants:
        depths = variant.allele_depths
        total = sum(depths) if depths else 0
        balance = (
            depths[-1] / total if depths and len(depths) >= 2 and total > 0 else None
        )
        depth = variant.depth
        if depth is None and depths:
            depth = float(total)
        rows.append(
            _one_hot(balance, BALANCE_BOUNDS)
            + _one_hot(depth, DEPTH_BOUNDS)
            + _one_hot(variant.genotype_quality, QUALITY_BOUNDS)
            + [float(variant.filter_status not in ("PASS", ".", ""))]
        )
    return np.array(rows, dtype=np.float64).reshape(len(rows), COLUMN_COUNT)


def compute_log_quality_ratio(genome: Genome) -> NDArray[np.float64]:
    """Fit quality log ratios using common calls as examples of true variants.

    Fit heterozygous and homozygous calls separately. A group without common variants
    receives a neutral log ratio of zero.
    """
    features = build_quality_features(genome.variants)
    is_common = np.array(
        [
            variant.allele_frequency is not None
            and variant.allele_frequency >= COMMON_FREQUENCY
            for variant in genome.variants
        ],
        dtype=bool,
    )
    result = np.zeros(len(genome.variants))
    for group in (genome.is_heterozygous, genome.is_homozygous):
        rows = np.flatnonzero(group)
        common = rows[is_common[rows]]
        if not len(common):
            continue
        ratio = fit_density_ratio(
            positive=features[common], negative=features[rows], regularisation=1.0
        )
        result[rows] = ratio.log_ratio(features[rows])
    return result
