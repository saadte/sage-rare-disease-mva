"""Fixed gene assignments for cross-validation."""

import csv
from collections.abc import Iterable, Mapping
from pathlib import Path

FOLD_COUNT = 5
DEFAULT_GENE_FOLDS = Path(__file__).resolve().parents[2] / "config" / "gene_folds.tsv"


def read_gene_folds(path: Path) -> dict[str, int]:
    """Read unique gene symbols and fold indices from a TSV with gene/fold columns."""
    folds = {}
    with path.open(encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            gene, fold = row["gene"], int(row["fold"])
            if not gene or not 0 <= fold < FOLD_COUNT or gene in folds:
                message = f"invalid or duplicate fold assignment for {gene!r}"
                raise ValueError(message)
            folds[gene] = fold
    if not folds:
        message = f"{path} has no gene fold assignments"
        raise ValueError(message)
    return folds


def assign_gene_folds(
    genes: Iterable[str], *, known_folds: Mapping[str, int] | None = None
) -> dict[str, int]:
    """Select fixed folds for the supplied symbols.

    Use the bundled assignments unless a reference table is supplied. Raise ValueError
    for unknown symbols; extend the reference table before using a new gene vocabulary.
    """
    folds = read_gene_folds(DEFAULT_GENE_FOLDS) if known_folds is None else known_folds
    requested = set(genes)
    missing = sorted(requested - folds.keys())
    if missing:
        message = f"gene fold assignments are missing for {', '.join(missing[:10])}"
        raise ValueError(message)
    return {gene: folds[gene] for gene in sorted(requested)}
