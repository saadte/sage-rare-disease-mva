import functools
import gzip
import re
from dataclasses import dataclass
from pathlib import Path

from mvapipe import __version__
from mvapipe.process import capture, run_pipeline
from mvapipe.references import References
from mvapipe.stage_runner import StageError, StageRunner
from mvapipe.variant_table import (
    CONSEQUENCE_COLUMNS,
    INFO_COLUMNS,
    OPTIONAL_FORMAT_COLUMNS,
    SITE_COLUMNS,
)

CANONICAL_CHROMOSOMES = [str(number) for number in range(1, 23)] + ["X", "Y", "MT"]
# Keep this cache release aligned with environment.yml and config/resources.env.
VEP_CACHE_VERSION = "116"
CLINVAR_TAGS = "INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNSIGCONF,INFO/GENEINFO"
MINIMUM_RECORDS_PER_CHUNK = 20_000
COMPRESSION_THREADS = "4"


@dataclass(frozen=True)
class Chunk:
    """Chromosomes annotated in one VEP run."""

    name: str
    chromosomes: tuple[str, ...]
    record_count: int


def choose_sample(vcf_path: Path, sample: str | None) -> str:
    """Select the requested sample, or the VCF's only sample if none is requested.

    Raise ValueError if the requested sample is absent or no unique sample can be chosen.
    """
    samples = capture(["bcftools", "query", "-l", str(vcf_path)]).split()
    if sample is None and len(samples) == 1:
        return samples[0]
    if sample is None or sample not in samples:
        message = f"choose one of the samples in {vcf_path} with --sample: {samples}"
        raise ValueError(message)
    return sample


def _count_records(vcf_path: Path) -> int:
    return int(capture(["bcftools", "index", "-n", str(vcf_path)]))


def index_vcf(vcf_path: Path, *, log_path: Path) -> None:
    run_pipeline([["bcftools", "index", "-f", "-t", str(vcf_path)]], log_path=log_path)


def _write_chromosome_names(names_path: Path) -> None:
    renames = [(f"chr{name}", name) for name in CANONICAL_CHROMOSOMES[:-1]]
    renames += [("chrM", "MT"), ("chrMT", "MT"), ("M", "MT")]
    names_path.write_text(
        "".join(f"{old}\t{new}\n" for old, new in renames), encoding="utf8"
    )


def normalise(
    log_path: Path,
    *,
    vcf_path: Path,
    sample: str,
    normalised_path: Path,
    references: References,
) -> None:
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
            # Fail on REF mismatches to avoid annotating against the wrong assembly.
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
    index_vcf(normalised_path, log_path=log_path)


def plan_chunks(record_counts: dict[str, int]) -> list[Chunk]:
    """Group chromosomes in genome order for VEP annotation.

    Close each chunk after reaching MINIMUM_RECORDS_PER_CHUNK records. Merge any final
    undersized group into the preceding chunk, or keep it alone if it is the only group.
    """
    groups: list[tuple[list[str], int]] = []
    chromosomes: list[str] = []
    record_count = 0
    for chromosome in CANONICAL_CHROMOSOMES:
        if chromosome not in record_counts:
            continue
        chromosomes.append(chromosome)
        record_count += record_counts[chromosome]
        if record_count >= MINIMUM_RECORDS_PER_CHUNK:
            groups.append((chromosomes, record_count))
            chromosomes, record_count = [], 0
    if chromosomes and groups:
        last_chromosomes, last_record_count = groups.pop()
        groups.append((last_chromosomes + chromosomes, last_record_count + record_count))
    elif chromosomes:
        groups.append((chromosomes, record_count))
    return [
        Chunk(
            name=group[0] if len(group) == 1 else f"{group[0]}_{group[-1]}",
            chromosomes=tuple(group),
            record_count=count,
        )
        for group, count in groups
    ]


def count_records_by_chromosome(vcf_path: Path) -> dict[str, int]:
    """Read nonzero chromosome record counts from the VCF index."""
    counts = {}
    for line in capture(["bcftools", "index", "-s", str(vcf_path)]).splitlines():
        chromosome, _, count_text = line.split("\t")
        count = int(count_text)
        if count > 0:
            counts[chromosome] = count
    return counts


def annotate_chunk(
    log_path: Path,
    *,
    normalised_path: Path,
    chunk: Chunk,
    chunk_path: Path,
    references: References,
    thread_count: int,
) -> None:
    """Annotate a chromosome chunk with VEP.

    Raise StageError if the output record count differs from the input count.
    """
    run_pipeline(
        [
            ["bcftools", "view", "-r", ",".join(chunk.chromosomes), str(normalised_path)],
            [
                "vep", "--offline", "--cache", "--dir_cache", str(references.vep_cache),
                "--cache_version", VEP_CACHE_VERSION, "--assembly", "GRCh38",
                "--species", "homo_sapiens", "--fasta", str(references.fasta),
                "--format", "vcf", "--vcf", "--compress_output", "bgzip",
                "--no_stats", "--force_overwrite",
                "--fork", str(thread_count), "--buffer_size", "20000",
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
    index_vcf(chunk_path, log_path=log_path)
    record_count = _count_records(chunk_path)
    if record_count != chunk.record_count:
        message = (
            f"chunk {chunk.name}: VEP wrote {record_count} "
            f"of {chunk.record_count} records"
        )
        raise StageError(message)


def merge_chunks(
    log_path: Path,
    *,
    chunk_paths: list[Path],
    annotated_path: Path,
    references: References,
) -> None:
    """Concatenate chunks and add ClinVar by exact allele match."""
    # Index the concatenated VCF before passing it to bcftools annotate.
    concatenated_path = annotated_path.with_name("concatenated.vcf.gz")
    run_pipeline(
        [[
            "bcftools", "concat", "-Oz", "-o", str(concatenated_path),
            "--threads", COMPRESSION_THREADS, *map(str, chunk_paths),
        ]],
        log_path=log_path,
    )  # fmt: skip
    index_vcf(concatenated_path, log_path=log_path)
    run_pipeline(
        [[
            "bcftools", "annotate", "-a", str(references.clinvar), "-c", CLINVAR_TAGS,
            "-Oz", "-o", str(annotated_path),
            "--threads", COMPRESSION_THREADS, str(concatenated_path),
        ]],
        log_path=log_path,
    )  # fmt: skip
    index_vcf(annotated_path, log_path=log_path)
    concatenated_path.unlink()
    concatenated_path.with_name(concatenated_path.name + ".tbi").unlink()


def extract_table(log_path: Path, *, annotated_path: Path, table_path: Path) -> None:
    """Extract one TSV row per variant and gene into a gzip-compressed table."""
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
    # Append rows as a second gzip member; readers concatenate it with the header.
    with table_path.open("ab") as table:
        run_pipeline(
            [
                ["bcftools", "+split-vep", str(annotated_path), "-d", "-f", row_format],
                ["gzip", "-c"],
            ],
            log_path=log_path,
            stdout=table,
        )


def build_variant_table(
    vcf_path: Path,
    *,
    sample: str | None,
    directory: Path,
    references: References,
    thread_count: int,
    runner: StageRunner,
) -> Path:
    """Normalise and annotate a VCF, returning the table path under `directory`.

    Reuse completed stages when their outputs exist and their recorded inputs match.
    Records include the sample, mvapipe version, VEP cache path and release, plus file
    paths, sizes and modification times for the VCF, FASTA, AlphaMissense scores,
    ClinVar VCF and AlphaMissense plugin.
    """
    sample = choose_sample(vcf_path, sample)
    normalised_path = directory / "normalised.vcf.gz"
    annotated_path = directory / "annotated.vcf.gz"
    table_path = directory / "variants.tsv.gz"
    chunk_directory = directory / "chunks"
    chunk_directory.mkdir(exist_ok=True)
    input_paths = [
        vcf_path,
        references.fasta,
        references.alphamissense,
        references.clinvar,
        references.vep_plugins / "AlphaMissense.pm",
    ]
    source = {
        "sample": sample,
        "pipeline_version": __version__,
        "vep_cache": [str(references.vep_cache.resolve()), VEP_CACHE_VERSION],
        "inputs": [describe_file(path) for path in input_paths],
    }

    runner.run(
        "normalise",
        functools.partial(
            normalise,
            vcf_path=vcf_path,
            sample=sample,
            normalised_path=normalised_path,
            references=references,
        ),
        outputs=_with_index(normalised_path),
        parameters=source,
    )
    chunks = plan_chunks(count_records_by_chromosome(normalised_path))
    chunk_paths = [chunk_directory / f"vep.{chunk.name}.vcf.gz" for chunk in chunks]
    for chunk, chunk_path in zip(chunks, chunk_paths, strict=True):
        runner.run(
            f"annotate_{chunk.name}",
            functools.partial(
                annotate_chunk,
                normalised_path=normalised_path,
                chunk=chunk,
                chunk_path=chunk_path,
                references=references,
                thread_count=thread_count,
            ),
            outputs=_with_index(chunk_path),
            parameters={**source, "record_count": chunk.record_count},
        )
    chunk_names = [chunk.name for chunk in chunks]
    runner.run(
        "merge",
        functools.partial(
            merge_chunks,
            chunk_paths=chunk_paths,
            annotated_path=annotated_path,
            references=references,
        ),
        outputs=_with_index(annotated_path),
        parameters={**source, "chunks": chunk_names},
    )
    runner.run(
        "extract",
        functools.partial(
            extract_table, annotated_path=annotated_path, table_path=table_path
        ),
        outputs=[table_path],
        parameters={**source, "chunks": chunk_names},
    )
    return table_path


def describe_file(path: Path) -> dict[str, object]:
    """Return the resolved path, size and modification time of a file."""
    status = path.stat()
    return {
        "path": str(path.resolve()),
        "bytes": status.st_size,
        "modified": status.st_mtime_ns,
    }


def _with_index(vcf_path: Path) -> tuple[Path, Path]:
    return vcf_path, vcf_path.with_name(vcf_path.name + ".tbi")
