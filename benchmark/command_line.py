import argparse
import dataclasses
import json
import os
import random
from collections.abc import Callable, Sequence
from datetime import UTC, datetime
from pathlib import Path

from benchmark import cases, evaluation, posterior
from benchmark.annotation import annotate_spike_variants
from benchmark.comparison import compare_top1, read_ranks
from mvapipe import __version__
from mvapipe.frequency import FrequencyMode
from mvapipe.gene_validity import load_gene_validity
from mvapipe.ontology import PhenotypeModel
from mvapipe.posterior.analysis import PosteriorConfig
from mvapipe.posterior.training import (
    SEED,
    load_training_data,
    read_annotated_alleles,
    read_published_cases,
    select_training_cases,
)
from mvapipe.references import References
from mvapipe.scoring_config import ScoringConfig
from mvapipe.variant_table import load_variants


class _Arguments(argparse.Namespace):
    handler: Callable[["_Arguments"], int]
    background: Path
    refs: Path
    out: Path
    genes_per_group: int
    seed: int
    cases: Path
    threads: int
    case_archive: Path | None
    case_table: Path | None
    training_table: Path | None
    gene_folds: Path | None
    frequency_mode: str
    phenotype: str
    visible_clinvar: bool
    use_quality: bool
    remove_causes: bool
    max_cases: int | None
    first: Path
    second: Path
    column: str
    equal_gene_weight: bool
    replicates: int


def _positive_integer(value: str) -> int:
    number = int(value)
    if number <= 0:
        message = "must be a positive integer"
        raise argparse.ArgumentTypeError(message)
    return number


def main(arguments: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="benchmark")
    commands = parser.add_subparsers(required=True, metavar="COMMAND")
    make = commands.add_parser("make-cases", help="draw benchmark cases")
    make.add_argument("--background", type=Path, required=True)
    make.add_argument("--refs", type=Path, required=True)
    make.add_argument("--out", type=Path, required=True)
    make.add_argument("--genes-per-group", type=_positive_integer, default=40)
    make.add_argument("--seed", type=int, required=True)
    make.set_defaults(handler=_make_cases)
    run = commands.add_parser("run", help="rank spike-in cases in a background genome")
    run.add_argument("--background", type=Path, required=True)
    run.add_argument("--cases", type=Path, required=True)
    run.add_argument("--refs", type=Path, required=True)
    run.add_argument("--out", type=Path, required=True)
    run.add_argument("--threads", type=_positive_integer, default=os.cpu_count())
    run.set_defaults(handler=_benchmark)

    fitted = commands.add_parser("posterior", help="evaluate held-out published cases")
    fitted.add_argument("--background", type=Path, required=True)
    fitted.add_argument("--refs", type=Path, required=True)
    fitted.add_argument("--out", type=Path, required=True)
    fitted.add_argument("--case-archive", type=Path)
    fitted.add_argument("--case-table", type=Path)
    fitted.add_argument("--training-table", type=Path)
    fitted.add_argument("--gene-folds", type=Path)
    fitted.add_argument(
        "--frequency-mode", choices=[m.value for m in FrequencyMode], default="all"
    )
    fitted.add_argument(
        "--phenotype",
        choices=[p.value for p in posterior.PhenotypeCondition],
        default="real",
    )
    fitted.add_argument("--visible-clinvar", action="store_true")
    fitted.add_argument("--use-quality", action="store_true")
    fitted.add_argument("--remove-causes", action="store_true")
    fitted.add_argument("--max-cases", type=_positive_integer)
    fitted.set_defaults(handler=_posterior_benchmark)

    compare = commands.add_parser("compare", help="compare paired ranks by causal gene")
    compare.add_argument("--first", type=Path, required=True)
    compare.add_argument("--second", type=Path, required=True)
    compare.add_argument("--column", choices=["rank", "gene_rank"], required=True)
    compare.add_argument("--equal-gene-weight", action="store_true")
    compare.add_argument("--replicates", type=_positive_integer, default=10_000)
    compare.add_argument("--seed", type=int, required=True)
    compare.add_argument("--out", type=Path, required=True)
    compare.set_defaults(handler=_compare)
    parsed = parser.parse_args(arguments, namespace=_Arguments())
    return parsed.handler(parsed)


def _make_cases(arguments: _Arguments) -> int:
    references = References(arguments.refs)
    background = load_variants(arguments.background)
    drawn = cases.make_cases(
        cases.load_clinvar_pool(references.clinvar),
        load_gene_validity(references.gencc),
        PhenotypeModel(
            obo_path=references.hpo_obo,
            genes_to_phenotype_path=references.genes_to_phenotype,
        ),
        {variant.allele for variant in background},
        genes_per_group=arguments.genes_per_group,
        random_generator=random.Random(arguments.seed),
    )
    cases.write_cases(drawn, arguments.out)
    print(f"{len(drawn)} cases written to {arguments.out}")
    return 0


def _benchmark(arguments: _Arguments) -> int:
    directory = arguments.out
    directory.mkdir(parents=True, exist_ok=True)
    references = References(arguments.refs)
    config = ScoringConfig()
    background = load_variants(arguments.background)
    drawn = cases.read_cases(arguments.cases)
    annotated = annotate_spike_variants(
        drawn, directory / "spikes", references, arguments.threads
    )
    results = evaluation.evaluate(
        drawn,
        annotated,
        background,
        load_gene_validity(references.gencc),
        PhenotypeModel(
            obo_path=references.hpo_obo,
            genes_to_phenotype_path=references.genes_to_phenotype,
        ),
        config,
    )
    summary = evaluation.summarise(results)
    evaluation.write_case_table(results, directory / "benchmark_cases.tsv")
    evaluated = {result.case.name for result in results}
    report = {
        "meta": {
            "pipeline_version": __version__,
            "finished": _now(),
            "background": str(arguments.background),
            "cases_file": str(arguments.cases),
            "cases": len(drawn),
            "cases_already_in_background": [
                case.name for case in drawn if case.name not in evaluated
            ],
            "scoring_config": dataclasses.asdict(config),
        },
        "summary": summary,
    }
    (directory / "benchmark_summary.json").write_text(
        json.dumps(report, indent=1), encoding="utf8"
    )
    markdown = evaluation.markdown_of_summary(summary)
    (directory / "benchmark_summary.md").write_text(markdown, encoding="utf-8")
    print(markdown)
    return 0


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _posterior_benchmark(arguments: _Arguments) -> int:
    references = References(arguments.refs)
    model = PhenotypeModel(
        obo_path=references.hpo_obo,
        genes_to_phenotype_path=references.genes_to_phenotype,
    )
    frequency_mode = FrequencyMode(arguments.frequency_mode)
    validity = load_gene_validity(references.gencc)
    training = load_training_data(
        references,
        validity,
        model,
        allele_table_path=arguments.training_table,
        gene_folds_path=arguments.gene_folds,
        frequency_mode=frequency_mode,
    )
    selected = training.cases
    if arguments.case_archive is not None:
        selected = select_training_cases(
            read_published_cases(arguments.case_archive), validity, model.ontology
        )
    if arguments.case_table is not None:
        extra = read_annotated_alleles(
            arguments.case_table, frequency_mode=frequency_mode
        )
        conflicting = [
            key
            for key in training.alleles.keys() & extra.keys()
            if training.alleles[key] != extra[key]
        ]
        if conflicting:
            message = f"training and case annotations disagree for {conflicting[0]}"
            raise ValueError(message)
        training = dataclasses.replace(
            training,
            alleles=training.alleles | extra,
        )
    selected = sorted(selected, key=lambda case: case.identifier)[: arguments.max_cases]
    if not selected:
        message = "no supported evaluation cases were selected"
        raise ValueError(message)
    options = posterior.EvaluationOptions(
        phenotype=posterior.PhenotypeCondition(arguments.phenotype),
        mask_clinvar=not arguments.visible_clinvar,
        remove_causes=arguments.remove_causes,
        posterior=PosteriorConfig(use_quality=arguments.use_quality),
    )
    results = posterior.evaluate_cases(
        selected,
        load_variants(arguments.background, frequency_mode=frequency_mode),
        training,
        options=options,
    )
    arguments.out.mkdir(parents=True, exist_ok=True)
    posterior.write_results(results, arguments.out / "cases.csv")
    report = {
        "pipeline_version": __version__,
        "selected_cases": len(results),
        "ranked_cases": sum(row.status == "ranked" for row in results),
        "phenotype": arguments.phenotype,
        "masked_clinvar": options.mask_clinvar,
        "remove_causes": options.remove_causes,
        "use_quality": arguments.use_quality,
        "frequency_mode": frequency_mode.value,
        "selection_seed": SEED,
        "background": str(arguments.background),
        "genotype_top1": sum(row.status == "ranked" and row.rank == 1 for row in results),
        "gene_top1": sum(
            row.status == "ranked" and row.gene_rank == 1 for row in results
        ),
    }
    (arguments.out / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(f"{len(results)} cases written to {arguments.out}")
    return 0


def _compare(arguments: _Arguments) -> int:
    result = compare_top1(
        read_ranks(arguments.first, column=arguments.column),
        read_ranks(arguments.second, column=arguments.column),
        equal_gene_weight=arguments.equal_gene_weight,
        replicates=arguments.replicates,
        seed=arguments.seed,
    )
    arguments.out.write_text(
        json.dumps(dataclasses.asdict(result), indent=2) + "\n", encoding="utf-8"
    )
    return 0
