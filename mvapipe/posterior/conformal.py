"""Candidate sets from held-out truth probabilities."""

import csv
import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from mvapipe.posterior.hypotheses import RankedAlleleSets
from mvapipe.variant import Allele


@dataclass(frozen=True)
class CandidateSet:
    """Allele sets above a calibrated cutoff, with the unexplained label kept separate."""

    allele_sets: tuple[frozenset[Allele], ...]
    includes_unexplained: bool
    minimum_probability: float
    alpha: float
    calibration_count: int


def read_truth_probabilities(path: Path) -> list[float]:
    """Read held-out truth probabilities from a CSV truth_probability column.

    Use records from the same method, phenotype and quality conditions as the target
    analysis. Include cause-removal controls when calibrating the unexplained label.
    """
    with path.open(encoding="utf-8", newline="") as source:
        probabilities = [
            float(row["truth_probability"]) for row in csv.DictReader(source)
        ]
    if not probabilities or any(
        not math.isfinite(value) or not 0 <= value <= 1 for value in probabilities
    ):
        message = f"{path} needs finite truth probabilities between 0 and 1"
        raise ValueError(message)
    return probabilities


def select_candidate_set(
    ranked: RankedAlleleSets, truth_probabilities: Sequence[float], *, alpha: float
) -> CandidateSet:
    """Select labels using the finite-sample split-conformal quantile.

    Supply held-out calibration probabilities and a ranking containing every allele
    set. Raise ValueError if truncation could hide a selected label. Coverage requires
    exchangeable calibration and test records; it is unverified for clinical use.
    """
    if not 0 < alpha < 1:
        message = "alpha must be strictly between 0 and 1"
        raise ValueError(message)
    if not truth_probabilities or any(
        not math.isfinite(value) or not 0 <= value <= 1 for value in truth_probabilities
    ):
        message = "calibration needs finite truth probabilities between 0 and 1"
        raise ValueError(message)
    ordered_scores = sorted(1 - value for value in truth_probabilities)
    count = len(ordered_scores)
    quantile_index = math.ceil((count + 1) * (1 - alpha))
    minimum = 0.0 if quantile_index > count else 1 - ordered_scores[quantile_index - 1]
    remaining = 1 - ranked.unexplained - math.fsum(ranked.posteriors.tolist())
    if not ranked.is_complete and remaining + 1e-12 >= minimum:
        message = "unreported probability mass could contain a candidate set member"
        raise ValueError(message)
    members = tuple(
        alleles
        for alleles, probability in zip(
            ranked.allele_sets, ranked.posteriors.tolist(), strict=True
        )
        if probability >= minimum
    )
    return CandidateSet(
        allele_sets=members,
        includes_unexplained=ranked.unexplained >= minimum,
        minimum_probability=minimum,
        alpha=alpha,
        calibration_count=count,
    )


def write_candidate_set(selection: CandidateSet, path: Path) -> None:
    """Write the selected allele sets and unexplained label as TSV rows."""
    with path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["label", "variant_1", "variant_2", "minimum_probability", "alpha", "n"]
        )
        calibration = [
            format(selection.minimum_probability, ".17g"),
            format(selection.alpha, ".17g"),
            selection.calibration_count,
        ]
        for alleles in selection.allele_sets:
            cells = [allele.to_string() for allele in sorted(alleles)]
            writer.writerow(
                [
                    "genotype",
                    *cells,
                    *([""] * (2 - len(cells))),
                    *calibration,
                ]
            )
        if selection.includes_unexplained:
            writer.writerow(["unexplained", "", "", *calibration])
