import argparse
import dataclasses
import json
import os
from collections.abc import Callable, Sequence
from datetime import UTC, datetime
from pathlib import Path

from mvapipe import __version__, submission
from mvapipe.gene_validity import load_gene_validity
from mvapipe.ontology import Phenotype_model, read_term_file
from mvapipe.ranked_table import write_ranked_table
from mvapipe.ranking import is_male, rank
from mvapipe.references import References
from mvapipe.scoring_config import Scoring_config
from mvapipe.stage_runner import Stage_runner
from mvapipe.variant_table import load_variants
from mvapipe.vcf_steps import build_variant_table

RANKED_TABLE_ROWS = 50


class _Arguments(argparse.Namespace):
    handler: Callable[["_Arguments"], int]
    out: Path
    refs: Path
    vcf: Path
    table: Path | None
    hpo: Path
    blind: bool
    sample: str | None
    proband_id: str
    label: str | None
    threads: int
    min_relative_score: float
    predictions: Path


def _positive_integer(value: str) -> int:
    number = int(value)
    if number <= 0:
        message = "must be a positive integer"
        raise argparse.ArgumentTypeError(message)
    return number


def _score_fraction(value: str) -> float:
    number = float(value)
    if not 0 <= number <= 1:
        message = "must be between 0 and 1"
        raise argparse.ArgumentTypeError(message)
    return number


def main(arguments: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="mvapipe")
    commands = parser.add_subparsers(required=True, metavar="COMMAND")

    run = commands.add_parser("run", help="rank a sample's variants; write predictions")
    source = run.add_mutually_exclusive_group(required=True)
    source.add_argument("--vcf", type=Path)
    source.add_argument("--table", type=Path, help="variants.tsv.gz of an earlier run")
    phenotype = run.add_mutually_exclusive_group(required=True)
    phenotype.add_argument("--hpo", type=Path, help="file with one HPO id per line")
    phenotype.add_argument(
        "--blind", action="store_true", help="rank without a phenotype"
    )
    run.add_argument("--out", type=Path, required=True, help="directory; a rerun resumes")
    run.add_argument("--refs", type=Path, required=True)
    run.add_argument("--sample", help="required when the VCF has several samples")
    run.add_argument("--proband-id", default="PROBAND01")
    run.add_argument(
        "--label", help="name of the output files (default: phenotype, blind)"
    )
    run.add_argument("--threads", type=_positive_integer, default=os.cpu_count())
    run.add_argument(
        "--min-relative-score",
        type=_score_fraction,
        default=0.0,
        help="after the best candidate, write only those scoring this fraction of it",
    )
    run.set_defaults(handler=_run)

    validate = commands.add_parser("validate", help="check a predictions file")
    validate.add_argument("predictions", type=Path)
    validate.add_argument("--proband-id", default="PROBAND01")
    validate.set_defaults(handler=_validate)

    parsed = parser.parse_args(arguments, namespace=_Arguments())
    return parsed.handler(parsed)


def _run(arguments: _Arguments) -> int:
    directory = arguments.out
    directory.mkdir(parents=True, exist_ok=True)
    references = References(arguments.refs)
    runner = Stage_runner(directory)
    table_path = arguments.table or build_variant_table(
        arguments.vcf,
        sample=arguments.sample,
        directory=directory,
        references=references,
        threads=arguments.threads,
        runner=runner,
    )
    variants, table_counts = load_variants(table_path)

    terms: list[str] = []
    phenotype_matches = None
    if not arguments.blind:
        model = Phenotype_model(references.hpo_obo, references.genes_to_phenotype)
        terms = read_term_file(arguments.hpo, model.ontology)
        genes = {variant.gene for variant in variants}
        phenotype_matches = {gene: model.gene_score(gene, terms) for gene in genes}
    config = Scoring_config()
    candidates = rank(
        variants, load_gene_validity(references.gencc), phenotype_matches, config
    )

    label = arguments.label or ("blind" if arguments.blind else "phenotype")
    predictions_path = directory / f"{label}.csv"
    write_ranked_table(candidates[:RANKED_TABLE_ROWS], directory / f"ranked.{label}.tsv")
    rows = submission.rows_of_candidates(
        candidates,
        proband_id=arguments.proband_id,
        minimum_relative_score=arguments.min_relative_score,
    )
    submission.write_predictions(rows, predictions_path)
    validation = submission.validate(predictions_path, proband_id=arguments.proband_id)
    manifest = {
        "pipeline_version": __version__,
        "finished": _now(),
        "table": str(table_path),
        "table_counts": table_counts,
        "stage_results": runner.results,
        "stage_seconds": runner.seconds,
        "phenotype_terms": len(terms),
        "is_male": is_male(variants),
        "candidates": len(candidates),
        "prediction_rows": len(rows),
        "minimum_relative_score": arguments.min_relative_score,
        "validation": dataclasses.asdict(validation),
        "scoring_config": dataclasses.asdict(config),
        "references": json.loads(references.manifest.read_text(encoding="utf8")),
    }
    (directory / f"manifest.{label}.json").write_text(
        json.dumps(manifest, indent=1), encoding="utf8"
    )
    print(f"{len(terms)} phenotype terms, {len(candidates)} candidates")
    print(f"{len(rows)} rows written to {predictions_path}")
    return _report(validation)


def _validate(arguments: _Arguments) -> int:
    validation = submission.validate(
        arguments.predictions, proband_id=arguments.proband_id
    )
    return _report(validation)


def _report(validation: submission.Validation) -> int:
    for error in validation.errors:
        print(f"error: {error}")
    for warning in validation.warnings:
        print(f"warning: {warning}")
    if validation.errors:
        return 1
    print("valid")
    return 0


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")
