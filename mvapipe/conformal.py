import math
from collections.abc import Sequence
from fractions import Fraction
from pathlib import Path


def read_scores(score_path: Path) -> list[float]:
    """Read finite, nonnegative calibration scores, one per nonblank line.

    Ignore text after '#'. Raise ValueError for an invalid score or a file
    containing no scores, with the line number included for invalid scores.
    """
    scores: list[float] = []
    with score_path.open(encoding="utf8") as lines:
        for number, raw_line in enumerate(lines, start=1):
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            try:
                score = float(line)
            except ValueError:
                message = f"{score_path}, line {number}: not a number"
                raise ValueError(message) from None
            if not math.isfinite(score) or score < 0:
                message = f"{score_path}, line {number}: not a score of 0 or more"
                raise ValueError(message)
            scores.append(score)
    if not scores:
        message = f"{score_path} has no scores"
        raise ValueError(message)
    return scores


def score_threshold(calibration_scores: Sequence[float], coverage: Fraction) -> float:
    """Return the split conformal threshold for the requested coverage.

    For n calibration scores, the threshold is the k-th highest score with
    k = ceil((n + 1) * coverage). Use 0 if k exceeds n, which admits all
    nonnegative candidate scores. The caller supplies finite, nonnegative
    scores and coverage between 0 and 1, excluding both endpoints.
    """
    ordered = sorted(calibration_scores)
    index = len(ordered) - math.ceil((len(ordered) + 1) * coverage)
    return ordered[index] if index >= 0 else 0.0
