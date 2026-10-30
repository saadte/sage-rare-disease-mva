import argparse
import csv
import gzip
from collections.abc import Sequence
from pathlib import Path

from mvapipe.variant_table import is_flagged_by_clinvar


def write_reduced_table(table_path: Path, reduced_path: Path) -> None:
    """Copy gene rows with non-MODIFIER impact and all ClinVar-flagged rows.

    Preserve the header and row order in the gzip-compressed output table.
    """
    with (
        gzip.open(table_path, "rt", encoding="utf8", newline="") as table,
        gzip.open(reduced_path, "wt", encoding="utf8", newline="") as reduced,
    ):
        rows = csv.DictReader(table, delimiter="\t")
        if rows.fieldnames is None:
            message = f"{table_path} has no header"
            raise ValueError(message)
        writer = csv.DictWriter(
            reduced, fieldnames=rows.fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            is_coding = row["SYMBOL"] not in (".", "") and row["IMPACT"] != "MODIFIER"
            if is_coding or is_flagged_by_clinvar(row):
                writer.writerow(row)


class _Arguments(argparse.Namespace):
    table: Path
    out: Path


def main(arguments: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="reduce-table")
    parser.add_argument("table", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parsed = parser.parse_args(arguments, namespace=_Arguments())
    write_reduced_table(parsed.table, parsed.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
