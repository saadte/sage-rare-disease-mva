"""Allele-matched frequencies from indexed GRCh38 population VCFs."""

import csv
import gzip
import re
from collections.abc import Iterable, Sequence
from pathlib import Path
from tempfile import TemporaryDirectory

from mvapipe.frequency import MATCHED_COLUMNS, maximum_frequency, parse_frequency
from mvapipe.process import capture, run_pipeline
from mvapipe.variant import Allele


def read_table_alleles(table_path: Path) -> set[Allele]:
    """Read unique normalised alleles from a compressed variant table."""
    with gzip.open(table_path, "rt", encoding="utf-8", newline="") as source:
        return {
            Allele(row["CHROM"], int(row["POS"]), row["REF"], row["ALT"])
            for row in csv.DictReader(source, delimiter="\t")
        }


def read_population_rows(
    lines: Iterable[str], alleles: set[Allele]
) -> dict[Allele, float | None]:
    """Match bcftools query rows by CHROM, POS, REF and the corresponding ALT.

    Each row contains these four fields and a Number=A frequency field. Preserve a
    reported zero. Reject ALT/frequency length mismatches and invalid fractions.
    """
    matched: dict[Allele, float | None] = {}
    for line in lines:
        chromosome, position, reference, alternates, frequencies = line.rstrip().split(
            "\t"
        )
        alts = alternates.split(",")
        values = frequencies.split(",")
        if frequencies == ".":
            values = ["."] * len(alts)
        for alternate, field in zip(alts, values, strict=True):
            value = parse_frequency(field)
            allele = Allele(
                chromosome.removeprefix("chr"), int(position), reference, alternate
            )
            if allele in alleles:
                matched[allele] = maximum_frequency((matched.get(allele), value))
    return matched


def _query_population(
    paths: Sequence[Path], alleles: set[Allele], directory: Path, *, log_path: Path
) -> dict[Allele, float | None]:
    matched: dict[Allele, float | None] = {}
    covered_chromosomes: set[str] = set()
    for path in paths:
        header = capture(["bcftools", "view", "--header-only", str(path)])
        if not re.search(r"##INFO=<ID=AF,Number=A,Type=Float[,>]", header):
            message = f"{path} needs a Number=A Float INFO/AF field"
            raise ValueError(message)
        contigs = re.findall(r"^##contig=<ID=([^,>]+)", header, flags=re.MULTILINE)
        names = {name.removeprefix("chr"): name for name in contigs}
        covered_chromosomes.update(names)
        regions = sorted(
            {(names[a.chromosome], a.position) for a in alleles if a.chromosome in names}
        )
        if not regions:
            continue
        regions_path = directory / "regions.tsv"
        regions_path.write_text(
            "".join(f"{chromosome}\t{position}\n" for chromosome, position in regions),
            encoding="utf-8",
        )
        queried = directory / "queried.tsv"
        with queried.open("wb") as destination:
            run_pipeline(
                [
                    [
                        "bcftools",
                        "query",
                        "--regions-file",
                        str(regions_path),
                        "--format",
                        "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n",
                        str(path),
                    ]
                ],
                log_path=log_path,
                stdout=destination,
            )
        with queried.open(encoding="utf-8") as source:
            for allele, value in read_population_rows(source, alleles).items():
                matched[allele] = maximum_frequency((matched.get(allele), value))
    missing = {allele.chromosome for allele in alleles} - covered_chromosomes
    if missing:
        message = f"population VCFs lack requested contigs: {', '.join(sorted(missing))}"
        raise ValueError(message)
    return matched


def prepare_population_frequencies(
    *,
    table_path: Path,
    exomes: Sequence[Path],
    genomes: Sequence[Path],
    output_path: Path,
) -> None:
    """Write a frequency lookup for every allele using supplied gnomAD VCF shards.

    Both sources must use GRCh38, left-normalised alleles and Number=A INFO/AF fractions.
    A missing match stays missing. These global AF values can differ from the maximum
    population-specific frequencies in a comparator database.
    """
    if not exomes or not genomes:
        message = "supply both exome and genome population VCFs"
        raise ValueError(message)
    if output_path.resolve() in {p.resolve() for p in [table_path, *exomes, *genomes]}:
        message = "write the frequency lookup to a different path"
        raise ValueError(message)
    alleles = read_table_alleles(table_path)
    if not alleles:
        message = f"{table_path} has no alleles"
        raise ValueError(message)
    with TemporaryDirectory(dir=output_path.parent) as temporary:
        directory = Path(temporary)
        log_path = output_path.with_suffix(".log")
        exome = _query_population(exomes, alleles, directory, log_path=log_path)
        genome = _query_population(genomes, alleles, directory, log_path=log_path)
        result = directory / "frequencies.tsv"
        with result.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
            writer.writerow(["CHROM", "POS", "REF", "ALT", *MATCHED_COLUMNS])
            for allele in sorted(alleles):
                values = (exome.get(allele), genome.get(allele))
                writer.writerow(
                    [
                        allele.chromosome,
                        allele.position,
                        allele.reference,
                        allele.alternate,
                        *(
                            "." if value is None else format(value, ".17g")
                            for value in values
                        ),
                    ]
                )
        result.replace(output_path)
