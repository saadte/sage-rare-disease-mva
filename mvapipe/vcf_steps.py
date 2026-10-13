import functools
import gzip
import json
import re
from dataclasses import dataclass
from pathlib import Path

from mvapipe import __version__
from mvapipe.process import capture, run, run_pipeline
from mvapipe.references import References
from mvapipe.stage_runner import Stage_error, Stage_result, Stage_runner
from mvapipe.variant_table import (
    CONSEQUENCE_COLUMNS,
    INFO_COLUMNS,
    OPTIONAL_FORMAT_COLUMNS,
    SITE_COLUMNS,
)

CANONICAL_CHROMOSOMES = [str(number) for number in range(1, 23)] + ["X", "Y", "MT"]
# Must match the VEP version in environment.yml and the cache in config/resources.env.
VEP_CACHE_VERSION = "116"
CLINVAR_TAGS = "INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNSIGCONF,INFO/GENEINFO"
MINIMUM_RECORDS_PER_CHUNK = 20_000
COMPRESSION_THREADS = "4"


@dataclass(frozen=True)
class Chunk:
    """Chromosomes annotated in one VEP run."""

    name: str
    chromosomes: tuple[str, ...]
    records: int


def sample_to_analyse(vcf_path: Path, sample: str | None) -> str:
    samples = capture(["bcftools", "query", "-l", str(vcf_path)]).split()
    if sample is None and len(samples) == 1:
        return samples[0]
    if sample is None or sample not in samples:
        message = f"choose one of the samples in {vcf_path} with --sample: {samples}"
        raise ValueError(message)
    return sample


def _record_count(vcf_path: Path) -> int:
    return int(capture(["bcftools", "index", "-n", str(vcf_path)]))


def _index(vcf_path: Path, log_path: Path) -> None:
    run(["bcftools", "index", "-f", "-t", str(vcf_path)], log_path=log_path)


def _write_chromosome_names(names_path: Path) -> None:
    renames = [(f"chr{name}", name) for name in CANONICAL_CHROMOSOMES[:-1]]
    renames += [("chrM", "MT"), ("chrMT", "MT"), ("M", "MT")]
    names_path.write_text(
        "".join(f"{old}\t{new}\n" for old, new in renames), encoding="utf8"
    )


def normalise(
    vcf_path: Path,
    sample: str,
    normalised_path: Path,
    references: References,
    log_path: Path,
) -> Stage_result:
    """Split and left-align one sample's alternate alleles on canonical chromosomes."""
    names_path = normalised_path.with_name("chromosome_names.txt")
    _write_chromosome_names(names_path)
    run_pipeline(
        [
            [
                "bcftools", "annotate", "--rename-chrs", str(names_path),
                "-x", "INFO", "-Ou", str(vcf_path),
            ],
            [
                "bcftools", "view", "-s", sample,
                "-t", ",".join(CANONICAL_CHROMOSOMES), "-Ou",
            ],
            # -c e: stop at the first REF allele that differs from the reference.
            [
                "bcftools", "norm", "-m", "-any", "-f", str(references.fasta), "-c", "e",
                "--old-rec-tag", "ORIG", "-Ou",
            ],
            [
                "bcftools", "view", "-i", 'GT="alt" && ALT!="*"', "-Oz",
                "-o", str(normalised_path), "--threads", COMPRESSION_THREADS,
            ],
        ],
        log_path=log_path,
    )  # fmt: skip
    _index(normalised_path, log_path)
    return {"records": _record_count(normalised_path)}


def plan_chunks(record_counts: dict[str, int]) -> list[Chunk]:
    """One chunk per chromosome; chromosomes below the minimum share a chunk."""
    groups: list[tuple[list[str], int]] = []
    chromosomes: list[str] = []
    records = 0
    for chromosome in CANONICAL_CHROMOSOMES:
        if chromosome not in record_counts:
            continue
        chromosomes.append(chromosome)
        records += record_counts[chromosome]
        if records >= MINIMUM_RECORDS_PER_CHUNK:
            groups.append((chromosomes, records))
            chromosomes, records = [], 0
    if chromosomes and groups:
        last_chromosomes, last_records = groups.pop()
        groups.append((last_chromosomes + chromosomes, last_records + records))
    elif chromosomes:
        groups.append((chromosomes, records))
    return [
        Chunk(
            name=group[0] if len(group) == 1 else f"{group[0]}_{group[-1]}",
            chromosomes=tuple(group),
            records=count,
        )
        for group, count in groups
    ]


def chromosome_record_counts(vcf_path: Path) -> dict[str, int]:
    counts = {}
    for line in capture(["bcftools", "index", "-s", str(vcf_path)]).splitlines():
        chromosome, _, records = line.split("\t")
        count = int(records)
        if count > 0:
            counts[chromosome] = count
    return counts


def annotate_chunk(
    normalised_path: Path,
    chunk: Chunk,
    chunk_path: Path,
    references: References,
    threads: int,
    log_path: Path,
) -> Stage_result:
    run_pipeline(
        [
            ["bcftools", "view", "-r", ",".join(chunk.chromosomes), str(normalised_path)],
            [
                "vep", "--offline", "--cache", "--dir_cache", str(references.vep_cache),
                "--cache_version", VEP_CACHE_VERSION, "--assembly", "GRCh38",
                "--species", "homo_sapiens", "--fasta", str(references.fasta),
                "--format", "vcf", "--vcf", "--compress_output", "bgzip",
                "--no_stats", "--force_overwrite",
                "--fork", str(threads), "--buffer_size", "20000",
                "--symbol", "--biotype", "--canonical", "--mane", "--numbers",
                "--pick_allele_gene",
                "--af_gnomade", "--af_gnomadg", "--max_af", "--check_existing",
                "--dir_plugins", str(references.vep_plugins),
                "--plugin", f"AlphaMissense,file={references.alphamissense}",
                "-o", str(chunk_path),
            ],
        ],
        log_path=log_path,
    )  # fmt: skip
    _index(chunk_path, log_path)
    records = _record_count(chunk_path)
    if records != chunk.records:
        message = f"chunk {chunk.name}: VEP wrote {records} of {chunk.records} records"
        raise Stage_error(message)
    return {"records": records}


def merge_chunks(
    chunk_paths: list[Path], annotated_path: Path, references: References, log_path: Path
) -> Stage_result:
    """Concatenate chunks and add ClinVar by exact allele match."""
    # bcftools annotate reads an indexed file.
    concatenated_path = annotated_path.with_name("concatenated.vcf.gz")
    run(
        [
            "bcftools", "concat", "-Oz", "-o", str(concatenated_path),
            "--threads", COMPRESSION_THREADS, *map(str, chunk_paths),
        ],
        log_path=log_path,
    )  # fmt: skip
    _index(concatenated_path, log_path)
    run(
        [
            "bcftools", "annotate", "-a", str(references.clinvar), "-c", CLINVAR_TAGS,
            "-Oz", "-o", str(annotated_path), "--threads", COMPRESSION_THREADS,
            str(concatenated_path),
        ],
        log_path=log_path,
    )  # fmt: skip
    _index(annotated_path, log_path)
    concatenated_path.unlink()
    concatenated_path.with_name(concatenated_path.name + ".tbi").unlink()
    return {"records": _record_count(annotated_path)}


def extract_table(annotated_path: Path, table_path: Path, log_path: Path) -> Stage_result:
    header = capture(["bcftools", "view", "-h", str(annotated_path)])
    format_tags = set(re.findall(r"##FORMAT=<ID=([^,>]+)", header))
    format_columns = ["GT"] + [
        tag for tag in OPTIONAL_FORMAT_COLUMNS if tag in format_tags
    ]
    columns = SITE_COLUMNS + INFO_COLUMNS + format_columns + CONSEQUENCE_COLUMNS
    row_format = (
        "\\t".join(f"%{column}" for column in SITE_COLUMNS)
        + "".join(f"\\t%INFO/{column}" for column in INFO_COLUMNS)
        + "["
        + "".join(f"\\t%{column}" for column in format_columns)
        + "]"
        + "".join(f"\\t%{column}" for column in CONSEQUENCE_COLUMNS)
        + "\\n"
    )
    with gzip.open(table_path, "wt", encoding="utf8", newline="") as table:
        table.write("\t".join(columns) + "\n")
    # The rows form a second gzip member, which readers join to the first.
    with table_path.open("ab") as table:
        run_pipeline(
            [
                ["bcftools", "+split-vep", str(annotated_path), "-d", "-f", row_format],
                ["gzip", "-c"],
            ],
            log_path=log_path,
            stdout=table,
        )
    return {}


def build_variant_table(
    vcf_path: Path,
    *,
    sample: str | None,
    directory: Path,
    references: References,
    threads: int,
    runner: Stage_runner,
) -> Path:
    sample = sample_to_analyse(vcf_path, sample)
    normalised_path = directory / "normalised.vcf.gz"
    annotated_path = directory / "annotated.vcf.gz"
    table_path = directory / "variants.tsv.gz"
    chunk_directory = directory / "chunks"
    chunk_directory.mkdir(exist_ok=True)
    status = vcf_path.stat()
    source = {
        "vcf": str(vcf_path.resolve()),
        "bytes": status.st_size,
        "modified": status.st_mtime_ns,
        "sample": sample,
        "pipeline_version": __version__,
        "references": json.loads(references.manifest.read_text(encoding="utf8")),
    }

    runner.run(
        "normalise",
        functools.partial(normalise, vcf_path, sample, normalised_path, references),
        outputs=_indexed_outputs(normalised_path),
        parameters=source,
    )
    chunks = plan_chunks(chromosome_record_counts(normalised_path))
    chunk_paths = [chunk_directory / f"vep.{chunk.name}.vcf.gz" for chunk in chunks]
    for chunk, chunk_path in zip(chunks, chunk_paths, strict=True):
        runner.run(
            f"annotate_{chunk.name}",
            functools.partial(
                annotate_chunk, normalised_path, chunk, chunk_path, references, threads
            ),
            outputs=_indexed_outputs(chunk_path),
            parameters={**source, "records": chunk.records},
        )
    chunk_names = [chunk.name for chunk in chunks]
    runner.run(
        "merge",
        functools.partial(merge_chunks, chunk_paths, annotated_path, references),
        outputs=_indexed_outputs(annotated_path),
        parameters={**source, "chunks": chunk_names},
    )
    runner.run(
        "extract",
        functools.partial(extract_table, annotated_path, table_path),
        outputs=[table_path],
        parameters={**source, "chunks": chunk_names},
    )
    return table_path


def _indexed_outputs(vcf_path: Path) -> tuple[Path, Path]:
    return vcf_path, vcf_path.with_name(vcf_path.name + ".tbi")
