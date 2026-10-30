import functools
from pathlib import Path

from mvapipe.posterior.training import read_published_cases, write_allele_vcf
from mvapipe.process import run_pipeline
from mvapipe.references import References
from mvapipe.stage_runner import StageRunner
from mvapipe.vcf_steps import build_variant_table, describe_file, index_vcf


def _write_checked_alleles(
    log_path: Path, *, references: References, checked_path: Path
) -> None:
    """Write published causal alleles, excluding REF mismatches with the reference."""
    vcf_path = checked_path.with_name("published_alleles.vcf")
    cases = read_published_cases(references.phenopackets)
    write_allele_vcf(
        (allele.key for case in cases for allele in case.small_alleles), vcf_path
    )
    run_pipeline(
        [
            [
                "bcftools", "norm", "-f", str(references.fasta), "-c", "x",
                "-Oz", "-o", str(checked_path), str(vcf_path),
            ],
        ],
        log_path=log_path,
    )  # fmt: skip
    index_vcf(checked_path, log_path=log_path)


def prepare_training_alleles(references: References, *, thread_count: int) -> Path:
    """Annotate published causal alleles and return the variant table path.

    Reuse completed stages when their inputs match, as in `build_variant_table`.
    """
    directory = references.training_directory
    directory.mkdir(parents=True, exist_ok=True)
    runner = StageRunner(directory)
    checked_path = directory / "published_alleles.checked.vcf.gz"
    runner.run(
        "check",
        functools.partial(
            _write_checked_alleles, references=references, checked_path=checked_path
        ),
        outputs=[checked_path, checked_path.with_name(checked_path.name + ".tbi")],
        parameters={
            "phenopackets": describe_file(references.phenopackets),
            "fasta": describe_file(references.fasta),
        },
    )
    return build_variant_table(
        checked_path,
        sample=None,
        directory=directory,
        references=references,
        thread_count=thread_count,
        runner=runner,
    )
