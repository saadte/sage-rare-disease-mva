"""Population allele frequencies expressed as fractions."""

import csv
import gzip
import math
import re
from collections.abc import Iterable, Mapping
from enum import Enum
from pathlib import Path
from tempfile import TemporaryDirectory

from mvapipe.variant import Allele

VEP_COLUMNS = ("gnomADe_AF", "gnomADg_AF", "MAX_AF")
MATCHED_COLUMNS = ("gnomAD_matched_exome_AF", "gnomAD_matched_genome_AF")


class FrequencyMode(Enum):
    """Population sources included when assigning an allele frequency."""

    ALL = "all"
    EXOME = "exome"


def maximum_frequency(values: Iterable[float | None]) -> float | None:
    """Return the highest reported fraction, preserving zero and missing values.

    Raise ValueError for nonfinite values or fractions outside [0, 1].
    """
    known = []
    for value in values:
        if value is None:
            continue
        if not math.isfinite(value) or not 0 <= value <= 1:
            message = f"allele frequency {value} is outside [0, 1]"
            raise ValueError(message)
        known.append(value)
    return max(known, default=None)


def parse_frequency(field: str) -> float | None:
    """Read the largest '&'- or ','-separated frequency, or None if missing."""
    return maximum_frequency(
        float(value) for value in re.split("[&,]", field) if value not in ("", ".")
    )


def frequency_of_row(row: Mapping[str, str], mode: FrequencyMode) -> float | None:
    """Read the highest frequency from the selected VEP and matched sources."""
    columns = (
        ("gnomADe_AF", MATCHED_COLUMNS[0])
        if mode is FrequencyMode.EXOME
        else VEP_COLUMNS + MATCHED_COLUMNS
    )
    return maximum_frequency(parse_frequency(row.get(column, ".")) for column in columns)


def add_matched_frequencies(
    *, table_path: Path, frequencies_path: Path, output_path: Path
) -> None:
    """Add allele-matched gnomAD columns to a gzip-compressed variant table.

    The frequency TSV requires CHROM, POS, REF, ALT and both MATCHED_COLUMNS. Include
    every table allele, using '.' for an unreported frequency. Preserve calls, VEP
    annotations and row order. Raise ValueError for missing or conflicting allele rows.
    """
    if output_path.resolve() in (table_path.resolve(), frequencies_path.resolve()):
        message = "write the annotated table to a different path"
        raise ValueError(message)
    frequencies: dict[Allele, tuple[float | None, float | None]] = {}
    with frequencies_path.open(encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        required = {"CHROM", "POS", "REF", "ALT", *MATCHED_COLUMNS}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            message = f"{frequencies_path} lacks required allele or frequency columns"
            raise ValueError(message)
        for row in reader:
            allele = Allele(
                row["CHROM"].removeprefix("chr"),
                int(row["POS"]),
                row["REF"].upper(),
                row["ALT"].upper(),
            )
            values = tuple(parse_frequency(row[column]) for column in MATCHED_COLUMNS)
            exome, genome = values
            if allele in frequencies and frequencies[allele] != (exome, genome):
                message = f"conflicting population frequencies for {allele.to_string()}"
                raise ValueError(message)
            frequencies[allele] = (exome, genome)
    with TemporaryDirectory(dir=output_path.parent) as directory:
        temporary = Path(directory) / "frequencies.tsv.gz"
        _write_matched_table(table_path, frequencies, temporary)
        temporary.replace(output_path)


def _write_matched_table(
    table_path: Path,
    frequencies: Mapping[Allele, tuple[float | None, float | None]],
    output_path: Path,
) -> None:
    with (
        gzip.open(table_path, "rt", encoding="utf-8", newline="") as source,
        gzip.open(output_path, "wt", encoding="utf-8", newline="") as destination,
    ):
        rows = csv.DictReader(source, delimiter="\t")
        if rows.fieldnames is None:
            message = f"{table_path} has no header"
            raise ValueError(message)
        columns = list(rows.fieldnames)
        columns.extend(column for column in MATCHED_COLUMNS if column not in columns)
        writer = csv.DictWriter(
            destination, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            allele = Allele(row["CHROM"], int(row["POS"]), row["REF"], row["ALT"])
            if allele not in frequencies:
                message = f"no frequency lookup row for {allele.to_string()}"
                raise ValueError(message)
            for column, value in zip(MATCHED_COLUMNS, frequencies[allele], strict=True):
                previous = parse_frequency(row.get(column, "."))
                combined = maximum_frequency((previous, value))
                row[column] = "." if combined is None else format(combined, ".17g")
            writer.writerow(row)
