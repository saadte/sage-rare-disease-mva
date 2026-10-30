"""Paired top-one comparisons with resampling by causal gene."""

import csv
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class RankRecord:
    case: str
    gene: str
    rank: int | None


@dataclass(frozen=True)
class PairedAccuracy:
    cases: int
    genes: int
    first_top1: int
    second_top1: int
    difference_pp: float
    interval_95_pp: tuple[float, float]
    weighting: str


def read_ranks(path: Path, *, column: str) -> list[RankRecord]:
    """Read case, gene and a selected rank column; an empty rank denotes a miss."""
    records = []
    with path.open(encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source):
            rank = int(row[column]) if row[column] else None
            if not row["case"] or not row["gene"] or (rank is not None and rank < 1):
                message = f"invalid case, gene or rank in {path}"
                raise ValueError(message)
            records.append(RankRecord(row["case"], row["gene"], rank))
    return records


def compare_top1(
    first: Sequence[RankRecord],
    second: Sequence[RankRecord],
    *,
    equal_gene_weight: bool = False,
    replicates: int = 10_000,
    seed: int,
) -> PairedAccuracy:
    """Compare identical case lists, using paired gene-cluster bootstrap intervals.

    Reject duplicate IDs, unmatched cases and inconsistent causal genes. Include misses
    in both denominators. Equal gene weighting averages within-gene accuracy differences.
    The interval describes variation among these genes on the supplied backgrounds.
    """
    if replicates < 2:
        message = "at least two bootstrap replicates are required"
        raise ValueError(message)
    first_by_id = {row.case: row for row in first}
    second_by_id = {row.case: row for row in second}
    if len(first_by_id) != len(first) or len(second_by_id) != len(second):
        message = "case IDs must be unique in each comparison table"
        raise ValueError(message)
    if not first_by_id or first_by_id.keys() != second_by_id.keys():
        message = "comparison tables must contain the same nonempty set of cases"
        raise ValueError(message)
    grouped: dict[str, list[int]] = defaultdict(list)
    for identifier, left in first_by_id.items():
        right = second_by_id[identifier]
        if left.gene != right.gene:
            message = f"causal genes disagree for {identifier}"
            raise ValueError(message)
        grouped[left.gene].append(int(left.rank == 1) - int(right.rank == 1))
    sizes = np.array([len(grouped[gene]) for gene in sorted(grouped)])
    differences = np.array([sum(grouped[gene]) for gene in sorted(grouped)])
    values = differences / sizes
    difference = (
        float(values.mean())
        if equal_gene_weight
        else float(differences.sum() / sizes.sum())
    )
    generator = np.random.default_rng(seed)
    sampled = np.empty(replicates)
    for index in range(replicates):
        indices = generator.integers(len(sizes), size=len(sizes))
        sampled[index] = (
            values[indices].mean()
            if equal_gene_weight
            else differences[indices].sum() / sizes[indices].sum()
        )
    lower, upper = np.quantile(sampled * 100, [0.025, 0.975])
    return PairedAccuracy(
        cases=len(first),
        genes=len(grouped),
        first_top1=sum(row.rank == 1 for row in first),
        second_top1=sum(row.rank == 1 for row in second),
        difference_pp=100 * difference,
        interval_95_pp=(float(lower), float(upper)),
        weighting="gene" if equal_gene_weight else "case",
    )
