from collections.abc import Sequence
from dataclasses import dataclass

from mvapipe.variant import Allele

RANK_POINTS = ((1, 100), (3, 50), (5, 25), (10, 10))
Prediction = tuple[frozenset[Allele], float]


@dataclass(frozen=True)
class Score:
    rank_points: float
    f_max: float
    full_match_rank: int | None
    partial_match_rank: int | None


def _points_of_rank(rank: int) -> int:
    return next((points for lowest_rank, points in RANK_POINTS if rank <= lowest_rank), 0)


def score_predictions(
    predictions: Sequence[Prediction], truth: frozenset[Allele]
) -> Score:
    """Compute challenge rank points and F-max from allele sets and EPCR values.

    Award half rank points for a partial match to a causal pair when no full match exists.
    Compute F-max over all distinct EPCR thresholds.
    """
    by_epcr = sorted(enumerate(predictions), key=lambda entry: (-entry[1][1], entry[0]))
    ranked = [alleles for _, (alleles, _) in by_epcr]
    full_match_rank = next(
        (rank for rank, alleles in enumerate(ranked, start=1) if alleles == truth), None
    )
    partial_match_rank = None
    is_pair = len(truth) == 2
    if is_pair and full_match_rank is None:
        partial_match_rank = next(
            (rank for rank, alleles in enumerate(ranked, start=1) if alleles & truth),
            None,
        )
    if full_match_rank is not None:
        rank_points = float(_points_of_rank(full_match_rank))
    elif partial_match_rank is not None:
        rank_points = 0.5 * _points_of_rank(partial_match_rank)
    else:
        rank_points = 0.0

    f_max = 0.0
    for threshold in sorted({epcr for _, epcr in predictions}, reverse=True):
        predicted: set[Allele] = set()
        for alleles, epcr in predictions:
            if epcr >= threshold:
                predicted |= alleles
        true_positives = len(predicted & truth)
        precision = true_positives / len(predicted) if predicted else 0.0
        recall = true_positives / len(truth) if truth else 0.0
        if precision + recall > 0:
            f_max = max(f_max, 2 * precision * recall / (precision + recall))
    return Score(rank_points, f_max, full_match_rank, partial_match_rank)
