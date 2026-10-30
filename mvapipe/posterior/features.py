from collections.abc import Sequence
from typing import assert_never

import numpy as np
from numpy.typing import NDArray

from mvapipe.ranking import (
    INFRAME,
    LOSS_OF_FUNCTION,
    SPLICE_REGION,
    TRUNCATING,
    ClinvarClass,
    classify_clinvar,
    is_in_last_exon,
    parse_review_stars,
)
from mvapipe.variant import Variant

CONSEQUENCES = (
    "truncating", "splice_site", "start_lost", "stop_lost", "inframe",
    "missense", "splice_region", "synonymous", "other",
)  # fmt: skip
# Upper bounds for positive-frequency bins; zero and missing values have separate bins.
FREQUENCY_SCHEMES = {
    "coarse": (1e-4, 1e-3, 1e-2),
    "default": (1e-5, 1e-4, 1e-3, 1e-2, 5e-2),
    "fine": (3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 5e-2),
}
CLINVAR_LEVELS = (
    "pathogenic_2", "pathogenic_1", "pathogenic_0", "conflicting_pathogenic",
    "uncertain", "benign_1", "benign_0", "absent",
)  # fmt: skip


def feature_names(frequency_bounds: Sequence[float]) -> tuple[str, ...]:
    """Return feature names in matrix column order."""
    return (
        *(f"consequence_{name}" for name in CONSEQUENCES),
        "truncating_last_exon",
        "missense_alphamissense", "missense_alphamissense_above_0.5",
        "missense_alphamissense_above_0.8", "missense_without_alphamissense",
        "frequency_zero", *(f"frequency_below_{bound:g}" for bound in frequency_bounds),
        "frequency_common", "frequency_missing",
        *(f"clinvar_{level}" for level in CLINVAR_LEVELS),
    )  # fmt: skip


def find_clinvar_columns(frequency_bounds: Sequence[float]) -> tuple[list[int], int]:
    """Return all ClinVar column indices and the index for ClinVar absence."""
    names = feature_names(frequency_bounds)
    columns = [index for index, name in enumerate(names) if name.startswith("clinvar_")]
    return columns, names.index("clinvar_absent")


def remove_clinvar(
    features: NDArray[np.float64], frequency_bounds: Sequence[float]
) -> NDArray[np.float64]:
    """Copy the feature matrix and mark every allele as absent from ClinVar."""
    columns, absent_column = find_clinvar_columns(frequency_bounds)
    unlisted = features.copy()
    unlisted[:, columns] = 0
    unlisted[:, absent_column] = 1
    return unlisted


def _classify_consequence(variant: Variant) -> str:
    terms = set(variant.consequence.split("&"))
    if terms & TRUNCATING:
        return "truncating"
    if terms & LOSS_OF_FUNCTION:
        return "splice_site"
    for name in ("start_lost", "stop_lost"):
        if name in terms:
            return name
    if terms & INFRAME:
        return "inframe"
    if "missense_variant" in terms:
        return "missense"
    if terms & SPLICE_REGION:
        return "splice_region"
    if "synonymous_variant" in terms:
        return "synonymous"
    return "other"


def _classify_clinvar_level(variant: Variant) -> str:
    stars = parse_review_stars(variant.clinvar_review_status)
    clinvar = classify_clinvar(
        significance=variant.clinvar_significance, conflict=variant.clinvar_conflict
    )
    match clinvar:
        case ClinvarClass.PATHOGENIC:
            return f"pathogenic_{min(stars, 2)}"
        case ClinvarClass.CONFLICTING_WITH_PATHOGENIC:
            return "conflicting_pathogenic"
        case ClinvarClass.BENIGN:
            return "benign_1" if stars >= 1 else "benign_0"
        case (
            ClinvarClass.CONFLICTING_MOSTLY_BENIGN | ClinvarClass.CONFLICTING_SOME_BENIGN
        ):
            return "benign_0"
        case ClinvarClass.CONFLICTING | ClinvarClass.UNCERTAIN:
            return "uncertain"
        case ClinvarClass.UNCLASSIFIED:
            return "absent"
        case _ as unreachable:
            assert_never(unreachable)


def _find_frequency_bin(frequency: float | None, bounds: Sequence[float]) -> int:
    """Return a frequency-bin index: zero, bounded positive bins, common or missing.

    Each bound belongs to the bin below it. Values above the last bound are common.
    """
    if frequency is None:
        return len(bounds) + 2
    if frequency == 0:
        return 0
    return 1 + int(np.searchsorted(bounds, frequency, side="left"))


def build_feature_matrix(
    variants: Sequence[Variant], frequency_bounds: Sequence[float]
) -> NDArray[np.float64]:
    """Build one row per variant, with columns ordered by `feature_names`."""
    names = feature_names(frequency_bounds)
    matrix = np.zeros((len(variants), len(names)))
    column = {name: index for index, name in enumerate(names)}
    frequency_columns = [
        column["frequency_zero"],
        *(column[f"frequency_below_{bound:g}"] for bound in frequency_bounds),
        column["frequency_common"],
        column["frequency_missing"],
    ]
    for row, variant in enumerate(variants):
        consequence = _classify_consequence(variant)
        matrix[row, column[f"consequence_{consequence}"]] = 1
        if consequence == "truncating" and is_in_last_exon(variant.exon):
            matrix[row, column["truncating_last_exon"]] = 1
        if consequence == "missense":
            score = variant.alphamissense
            if score is None:
                matrix[row, column["missense_without_alphamissense"]] = 1
            else:
                matrix[row, column["missense_alphamissense"]] = score
                matrix[row, column["missense_alphamissense_above_0.5"]] = max(
                    score - 0.5, 0
                )
                matrix[row, column["missense_alphamissense_above_0.8"]] = max(
                    score - 0.8, 0
                )
        frequency_bin = _find_frequency_bin(variant.allele_frequency, frequency_bounds)
        matrix[row, frequency_columns[frequency_bin]] = 1
        matrix[row, column[f"clinvar_{_classify_clinvar_level(variant)}"]] = 1
    return matrix
