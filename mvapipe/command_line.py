import argparse
import os
from collections.abc import Callable, Sequence
from pathlib import Path
from tempfile import TemporaryDirectory

from mvapipe import submission
from mvapipe.frequency import FrequencyMode, add_matched_frequencies
from mvapipe.gene_validity import GeneValidity, load_gene_validity
from mvapipe.ontology import PhenotypeModel, read_term_file
from mvapipe.population import prepare_population_frequencies
from mvapipe.posterior.analysis import (
    PosteriorConfig,
    analyse_genome,
    convert_allele_sets,
    write_posterior_table,
)
from mvapipe.posterior.conformal import (
    read_truth_probabilities,
    select_candidate_set,
    write_candidate_set,
)
from mvapipe.posterior.preparation import prepare_training_alleles
from mvapipe.posterior.training import load_training_data
from mvapipe.ranked_table import write_ranked_table
from mvapipe.ranking import rank
from mvapipe.references import References
from mvapipe.scoring_config import ScoringConfig
from mvapipe.stage_runner import StageRunner
from mvapipe.variant import Variant
from mvapipe.variant_table import load_variants
from mvapipe.vcf_steps import build_variant_table

RANKED_TABLE_ROWS = 50
METHODS = ("score", "posterior")


class _Arguments(argparse.Namespace):
    handler: Callable[["_Arguments"], int]
    out: Path
    refs: Path
    vcf: Path
    table: Path | None
    hpo: Path
    blind: bool
    method: str
    sample: str | None
    proband_id: str
    label: str | None
    thread_count: int
    min_relative_score: float
    predictions: Path
    frequency_mode: str
    use_quality: bool
    unexplained_prior: float
    unlisted_prior: float
    calibration_probabilities: Path | None
    alpha: float
    gene_folds: Path | None
    training_table: Path | None
    source_table: Path
    frequencies: Path
    exomes: list[Path]
    genomes: list[Path]


def _parse_positive_integer(value: str) -> int:
    number = int(value)
    if number <= 0:
        message = "must be a positive integer"
        raise argparse.ArgumentTypeError(message)
    return number


def _parse_score_fraction(value: str) -> float:
    number = float(value)
    if not 0 <= number <= 1:
        message = "must be between 0 and 1"
        raise argparse.ArgumentTypeError(message)
    return number


def _parse_probability(value: str) -> float:
    number = float(value)
    if not 0 < number < 1:
        message = "must be strictly between 0 and 1"
        raise argparse.ArgumentTypeError(message)
    return number


def _add_thread_argument(parser: argparse.ArgumentParser) -> None:
    cpu_count = os.cpu_count()
    parser.add_argument(
        "--threads",
        dest="thread_count",
        type=_parse_positive_integer,
        default=cpu_count,
        required=cpu_count is None,
    )


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
    run.add_argument(
        "--method",
        choices=METHODS,
        default="score",
        help="rank by the score or by the posterior probability (default: score)",
    )
    run.add_argument("--out", type=Path, required=True, help="directory; a rerun resumes")
    run.add_argument("--refs", type=Path, required=True)
    run.add_argument(
        "--frequency-mode",
        choices=[mode.value for mode in FrequencyMode],
        default="all",
        help="use all population frequencies or exome frequencies only",
    )
    run.add_argument(
        "--use-quality", action="store_true", help="fit posterior call-quality evidence"
    )
    run.add_argument("--unexplained-prior", type=_parse_probability, default=0.5)
    run.add_argument("--unlisted-prior", type=_parse_probability, default=0.75)
    run.add_argument(
        "--training-table",
        type=Path,
        help="annotated causal alleles with matched frequencies",
    )
    run.add_argument("--gene-folds", type=Path, help="TSV of fixed gene/fold assignments")
    run.add_argument(
        "--calibration-probabilities",
        type=Path,
        help="CSV of held-out truth probabilities for a conformal candidate set",
    )
    run.add_argument("--alpha", type=_parse_probability, default=0.1)
    run.add_argument("--sample", help="required when the VCF has several samples")
    run.add_argument("--proband-id", default="PROBAND01")
    run.add_argument(
        "--label", help="name of the output files (default: phenotype, blind)"
    )
    _add_thread_argument(run)
    run.add_argument(
        "--min-relative-score",
        type=_parse_score_fraction,
        default=0.0,
        help="after the best candidate, write only those scoring this fraction of it",
    )
    run.set_defaults(handler=_run)

    prepare = commands.add_parser(
        "prepare-training",
        help="annotate the published causal alleles used by --method posterior",
    )
    prepare.add_argument("--refs", type=Path, required=True)
    _add_thread_argument(prepare)
    prepare.set_defaults(handler=_prepare_training)

    frequencies = commands.add_parser(
        "add-frequencies", help="add matched population frequencies to a variant table"
    )
    frequencies.add_argument("--table", dest="source_table", type=Path, required=True)
    frequencies.add_argument("--frequencies", type=Path, required=True)
    frequencies.add_argument("--out", type=Path, required=True)
    frequencies.set_defaults(handler=_add_frequencies)

    lookup = commands.add_parser(
        "prepare-frequencies", help="match alleles to indexed GRCh38 gnomAD VCFs"
    )
    lookup.add_argument("--table", dest="source_table", type=Path, required=True)
    lookup.add_argument("--exomes", type=Path, nargs="+", required=True)
    lookup.add_argument("--genomes", type=Path, nargs="+", required=True)
    lookup.add_argument("--out", type=Path, required=True)
    lookup.set_defaults(handler=_prepare_frequencies)

    validate = commands.add_parser("validate", help="check a predictions file")
    validate.add_argument("predictions", type=Path)
    validate.add_argument("--proband-id", default="PROBAND01")
    validate.set_defaults(handler=_validate)

    parsed = parser.parse_args(arguments, namespace=_Arguments())
    return parsed.handler(parsed)


def _run(arguments: _Arguments) -> int:
    if arguments.method != "posterior" and (
        arguments.use_quality
        or arguments.calibration_probabilities is not None
        or arguments.gene_folds is not None
        or arguments.training_table is not None
        or arguments.unexplained_prior != 0.5
        or arguments.unlisted_prior != 0.75
        or arguments.alpha != 0.1
    ):
        message = (
            "quality, training, calibration and prior options require --method posterior"
        )
        raise ValueError(message)
    if arguments.calibration_probabilities is not None and arguments.min_relative_score:
        message = "relative-score filtering would discard calibrated set members"
        raise ValueError(message)
    if arguments.calibration_probabilities is None and arguments.alpha != 0.1:
        message = "--alpha requires --calibration-probabilities"
        raise ValueError(message)
    directory = arguments.out
    directory.mkdir(parents=True, exist_ok=True)
    references = References(arguments.refs)
    table_path = arguments.table or build_variant_table(
        arguments.vcf,
        sample=arguments.sample,
        directory=directory,
        references=references,
        thread_count=arguments.thread_count,
        runner=StageRunner(directory),
    )
    variants = load_variants(
        table_path, frequency_mode=FrequencyMode(arguments.frequency_mode)
    )
    gene_validity = load_gene_validity(references.gencc)

    label = arguments.label or ("blind" if arguments.blind else "phenotype")
    if arguments.method == "score":
        predictions = _predict_by_score(
            arguments,
            references,
            variants,
            gene_validity,
            table_path=directory / f"ranked.{label}.tsv",
        )
    else:
        predictions = _predict_by_posterior(
            arguments,
            references,
            variants,
            gene_validity,
            table_path=directory / f"posterior.{label}.tsv",
        )
    predictions_path = directory / f"{label}.csv"
    rows = submission.build_prediction_rows(
        predictions,
        proband_id=arguments.proband_id,
        minimum_relative_score=arguments.min_relative_score,
    )
    with TemporaryDirectory(dir=directory) as temporary_directory:
        temporary = Path(temporary_directory) / predictions_path.name
        submission.write_predictions(rows, temporary)
        validation = submission.validate(temporary, proband_id=arguments.proband_id)
        if validation.errors:
            return _report(validation)
        temporary.replace(predictions_path)
    print(f"{len(rows)} rows written to {predictions_path}")
    return _report(validation)


def _load_phenotype_model(references: References) -> PhenotypeModel:
    return PhenotypeModel(
        obo_path=references.hpo_obo,
        genes_to_phenotype_path=references.genes_to_phenotype,
    )


def _predict_by_score(
    arguments: _Arguments,
    references: References,
    variants: Sequence[Variant],
    gene_validity: dict[str, GeneValidity],
    *,
    table_path: Path,
) -> list[submission.Prediction]:
    phenotype_matches = None
    if not arguments.blind:
        model = _load_phenotype_model(references)
        terms = read_term_file(arguments.hpo, model.ontology)
        print(f"{len(terms)} phenotype terms")
        genes = {variant.gene for variant in variants}
        phenotype_matches = {gene: model.score_gene(gene, terms) for gene in genes}
    candidates = rank(
        variants,
        gene_validity,
        phenotype_matches=phenotype_matches,
        config=ScoringConfig(),
    )
    write_ranked_table(candidates[:RANKED_TABLE_ROWS], table_path)
    return submission.convert_candidates(candidates)


def _predict_by_posterior(
    arguments: _Arguments,
    references: References,
    variants: Sequence[Variant],
    gene_validity: dict[str, GeneValidity],
    *,
    table_path: Path,
) -> list[submission.Prediction]:
    model = _load_phenotype_model(references)
    terms = None
    if not arguments.blind:
        terms = read_term_file(arguments.hpo, model.ontology)
        print(f"{len(terms)} phenotype terms")
    training = load_training_data(
        references,
        gene_validity,
        model,
        allele_table_path=arguments.training_table,
        gene_folds_path=arguments.gene_folds,
        frequency_mode=FrequencyMode(arguments.frequency_mode),
    )
    analysis = analyse_genome(
        variants,
        terms,
        training,
        top=None
        if arguments.calibration_probabilities is not None
        else RANKED_TABLE_ROWS,
        config=PosteriorConfig(
            unexplained_prior=arguments.unexplained_prior,
            unlisted_prior=arguments.unlisted_prior,
            use_quality=arguments.use_quality,
        ),
    )
    write_posterior_table(analysis, table_path, limit=RANKED_TABLE_ROWS)
    predictions = convert_allele_sets(analysis)
    if arguments.calibration_probabilities is None:
        return predictions
    selection = select_candidate_set(
        analysis.ranked,
        read_truth_probabilities(arguments.calibration_probabilities),
        alpha=arguments.alpha,
    )
    label = arguments.label or ("blind" if arguments.blind else "phenotype")
    write_candidate_set(selection, arguments.out / f"set.{label}.tsv")
    selected = set(selection.allele_sets)
    predictions = [
        prediction for prediction in predictions if prediction.alleles in selected
    ]
    if not 1 <= len(predictions) <= submission.MAX_ROWS:
        message = "the candidate set must contain 1 to 10 genotypes for Track 1 export"
        raise ValueError(message)
    return predictions


def _add_frequencies(arguments: _Arguments) -> int:
    add_matched_frequencies(
        table_path=arguments.source_table,
        frequencies_path=arguments.frequencies,
        output_path=arguments.out,
    )
    return 0


def _prepare_frequencies(arguments: _Arguments) -> int:
    prepare_population_frequencies(
        table_path=arguments.source_table,
        exomes=arguments.exomes,
        genomes=arguments.genomes,
        output_path=arguments.out,
    )
    return 0


def _prepare_training(arguments: _Arguments) -> int:
    table_path = prepare_training_alleles(
        References(arguments.refs), thread_count=arguments.thread_count
    )
    print(f"published alleles annotated in {table_path}")
    return 0


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
