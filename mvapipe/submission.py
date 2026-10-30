import csv
import io
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from mvapipe.ranking import Candidate
from mvapipe.variant import Allele

COLUMNS = [
    "proband_id",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "epcr",
    "finding_type",
    "notes",
]
REQUIRED_COLUMNS = COLUMNS[:10]
MAX_ROWS = 10
BYTE_ORDER_MARK = b"\xef\xbb\xbf"

Row = dict[str, str]


@dataclass
class Validation:
    """Validation errors and warnings; any error makes the predictions file invalid."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def format_chromosome(chromosome: str) -> str:
    """Add the submission's 'chr' prefix and map mitochondrial names to 'chrM'."""
    name = chromosome.removeprefix("chr")
    return "chrM" if name in ("MT", "M") else f"chr{name}"


@dataclass(frozen=True)
class Prediction:
    """A candidate genotype with scores and a note for the predictions file."""

    alleles: frozenset[Allele]
    score: float
    epcr: float
    note: str


def convert_candidates(candidates: Sequence[Candidate]) -> list[Prediction]:
    """Convert ranked candidates to predictions with gene and model notes."""
    return [
        Prediction(
            alleles=candidate.hypothesis.alleles,
            score=candidate.score,
            epcr=candidate.epcr,
            note=f"{candidate.hypothesis.gene} {candidate.hypothesis.model.value}",
        )
        for candidate in candidates
    ]


def build_prediction_rows(
    predictions: Sequence[Prediction], *, proband_id: str, minimum_relative_score: float
) -> list[Row]:
    """Build submission rows from predictions supplied in descending score order.

    Keep the first prediction and up to MAX_ROWS rows in total. Stop when a later score
    falls below `minimum_relative_score` times the first score.
    """
    rows: list[Row] = []
    for prediction in predictions[:MAX_ROWS]:
        if rows and prediction.score < minimum_relative_score * predictions[0].score:
            break
        row = dict.fromkeys(COLUMNS, "")
        row["proband_id"] = proband_id
        row["epcr"] = format(prediction.epcr, ".17g")
        row["finding_type"] = "primary"
        row["notes"] = prediction.note
        for index, allele in enumerate(sorted(prediction.alleles), start=1):
            row[f"chrom_{index}"] = format_chromosome(allele.chromosome)
            row[f"pos_{index}"] = str(allele.position)
            row[f"ref_{index}"] = allele.reference.upper()
            row[f"alt_{index}"] = allele.alternate.upper()
        rows.append(row)
    return rows


def write_predictions(rows: Sequence[Row], predictions_path: Path) -> None:
    with predictions_path.open("w", encoding="utf-8", newline="") as predictions:
        writer = csv.DictWriter(predictions, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _parse_allele(row: Row, index: int) -> Allele | None:
    chromosome, position, reference, alternate = (
        row[f"{field}_{index}"].strip() for field in ("chrom", "pos", "ref", "alt")
    )
    if not any((chromosome, position, reference, alternate)):
        return None
    return Allele(chromosome, int(position), reference.upper(), alternate.upper())


def _parse_finding_type(row: Row) -> str:
    return row.get("finding_type", "").strip().lower() or "primary"


def _find_allele_errors(alleles: Sequence[Allele]) -> list[str]:
    errors = []
    for allele in alleles:
        if not allele.chromosome.startswith("chr"):
            errors.append(f"chromosome {allele.chromosome!r} lacks the 'chr' prefix")
        if allele.position <= 0:
            errors.append("position must be a positive integer")
        bases = allele.reference + allele.alternate
        if not allele.reference or not allele.alternate or set(bases) - set("ACGTN"):
            errors.append(
                f"REF/ALT must be non-empty A/C/G/T/N strings, "
                f"got {allele.reference!r}>{allele.alternate!r}"
            )
    return errors


def _check_epcr_across_rows(
    rows: Sequence[Row], epcr_values: Sequence[float], validation: Validation
) -> None:
    if epcr_values.count(max(epcr_values)) > 1:
        validation.errors.append(
            "the highest epcr is shared by more than one row, which lowers F-max"
        )
    if list(epcr_values) != sorted(epcr_values, reverse=True):
        validation.warnings.append("rows are not sorted by epcr, highest first")
    finding_types = [_parse_finding_type(row) for row in rows]
    if "secondary" not in finding_types:
        return
    first_secondary = epcr_values[finding_types.index("secondary")]
    primary_values = [
        epcr
        for epcr, finding_type in zip(epcr_values, finding_types, strict=True)
        if finding_type == "primary"
    ]
    if first_secondary >= max(primary_values, default=0.0):
        validation.errors.append(
            "a secondary finding has an epcr at or above the top primary row; "
            "the scorer ignores finding_type"
        )


def _parse_epcr(row: Row, line_number: int, errors: list[str]) -> float | None:
    try:
        epcr = float(row["epcr"])
    except ValueError:
        errors.append(f"line {line_number}: epcr is not a number")
        return None
    if not 0 < epcr <= 1:
        errors.append(f"line {line_number}: epcr {row['epcr']} is outside (0, 1]")
    return epcr


def _check_rows(rows: Sequence[Row], validation: Validation) -> None:
    errors, warnings = validation.errors, validation.warnings
    epcr_values: list[float] = []
    seen: set[frozenset[Allele]] = set()
    for line_number, row in enumerate(rows, start=2):
        try:
            first, second = _parse_allele(row, 1), _parse_allele(row, 2)
        except ValueError:
            errors.append(f"line {line_number}: position is not an integer")
            continue
        if first is None:
            errors.append(f"line {line_number}: first variant is missing")
            continue
        alleles = [first] if second is None else [first, second]
        errors.extend(
            f"line {line_number}: {error}" for error in _find_allele_errors(alleles)
        )
        if first == second:
            warnings.append(f"line {line_number}: both variants in the row are identical")
        if frozenset(alleles) in seen:
            warnings.append(f"line {line_number}: duplicate of an earlier row")
        seen.add(frozenset(alleles))
        epcr = _parse_epcr(row, line_number, errors)
        if epcr is not None:
            epcr_values.append(epcr)
        if _parse_finding_type(row) not in ("primary", "secondary"):
            errors.append(
                f"line {line_number}: finding_type must be primary or secondary"
            )
    if len(epcr_values) == len(rows):
        _check_epcr_across_rows(rows, epcr_values, validation)


def _read_rows(text: str, errors: list[str]) -> list[Row]:
    reader = csv.reader(io.StringIO(text, newline=""))
    header = next(reader, [])
    missing = [column for column in REQUIRED_COLUMNS if column not in header]
    if missing:
        errors.append(f"missing column(s): {', '.join(missing)}")
        return []
    if len(header) != len(set(header)):
        errors.append("duplicate column names")
        return []
    rows: list[Row] = []
    for line_number, fields in enumerate(reader, start=2):
        if not fields:
            continue
        if len(fields) != len(header):
            errors.append(
                f"line {line_number}: expected {len(header)} fields, got {len(fields)}"
            )
        else:
            rows.append(dict(zip(header, fields, strict=True)))
    return rows


def validate(predictions_path: Path, *, proband_id: str) -> Validation:
    """Validate a predictions CSV against the Track 1 submission rules."""
    validation = Validation()
    errors = validation.errors
    content = predictions_path.read_bytes()
    if content.startswith(BYTE_ORDER_MARK):
        errors.append(
            "file starts with a UTF-8 byte-order mark, which the scorer rejects"
        )
    try:
        text = content.decode("utf-8-sig")
    except UnicodeDecodeError:
        errors.append("file is not valid UTF-8")
        return validation
    previous_error_count = len(errors)
    rows = _read_rows(text, errors)
    if len(errors) > previous_error_count:
        return validation
    if not rows:
        errors.append("no data rows")
        return validation
    if len(rows) > MAX_ROWS:
        errors.append(f"{len(rows)} rows; the maximum is {MAX_ROWS}")
    proband_ids = {row["proband_id"].strip() for row in rows}
    if proband_ids != {proband_id}:
        errors.append(
            f"proband_id must be {proband_id!r} on every row, found {sorted(proband_ids)}"
        )
    _check_rows(rows, validation)
    return validation
