import csv
import itertools
from collections import defaultdict
from collections.abc import Callable, Collection, Mapping, Sequence
from dataclasses import dataclass, replace
from enum import Enum
from pathlib import Path

from benchmark.cases import Case, Genotype
from benchmark.scoring import score_predictions
from mvapipe import submission
from mvapipe.gene_validity import GeneValidity
from mvapipe.ontology import PhenotypeMatch, PhenotypeModel
from mvapipe.ranking import (
    Candidate,
    Evidence,
    build_hypotheses,
    is_male,
    rank_hypotheses,
)
from mvapipe.scoring_config import ScoringConfig
from mvapipe.variant import Allele, Variant

SUMMARY_COLUMNS = [
    "clinvar",
    "phenotype",
    "cases",
    "n",
    "top1_pct",
    "top3_pct",
    "top5_pct",
    "top10_pct",
    "mean_rank_points",
    "mean_f_max",
]


class ClinvarCondition(Enum):
    VISIBLE = "visible"
    # Clear ClinVar fields on spiked variants.
    MASKED = "masked"


class PhenotypeCondition(Enum):
    EXACT = "exact"
    NOISY = "noisy"
    BLIND = "blind"


def _terms_of_case(case: Case, condition: PhenotypeCondition) -> tuple[str, ...] | None:
    match condition:
        case PhenotypeCondition.EXACT:
            return case.exact_terms
        case PhenotypeCondition.NOISY:
            return case.noisy_terms
        case PhenotypeCondition.BLIND:
            return None


@dataclass(frozen=True)
class Outcome:
    # Rank of the exact spiked allele set, or None if it is absent.
    rank: int | None
    gene_rank: int | None
    # Challenge scores computed from the top ten candidates.
    rank_points: float
    f_max: float


NOT_RANKED = Outcome(rank=None, gene_rank=None, rank_points=0.0, f_max=0.0)


@dataclass(frozen=True)
class CaseResult:
    case: Case
    outcomes: dict[tuple[ClinvarCondition, PhenotypeCondition], Outcome]


CASE_GROUPS: tuple[tuple[str, Callable[[Case], bool]], ...] = (
    ("all", lambda _: True),
    ("comphet", lambda case: case.genotype is Genotype.COMPOUND_HETEROZYGOUS),
    (
        "comphet lof+missense",
        lambda case: (
            case.genotype is Genotype.COMPOUND_HETEROZYGOUS
            and case.kinds == "lof+missense"
        ),
    ),
    ("hom", lambda case: case.genotype is Genotype.HOMOZYGOUS),
    ("het (dominant)", lambda case: case.genotype is Genotype.HETEROZYGOUS),
)


def _is_loaded_without_clinvar(variant: Variant) -> bool:
    return variant.vep_impact != "MODIFIER" and not variant.is_gene_from_clinvar


def spiked_variants(
    case: Case, annotated: Mapping[str, Sequence[Variant]], condition: ClinvarCondition
) -> list[Variant] | None:
    """Build annotated spike variants with synthetic genotypes and call quality.

    In masked mode, clear ClinVar fields and drop annotations that depend on them.
    Return None if any spike allele has no usable annotation.
    """
    is_homozygous = case.genotype is Genotype.HOMOZYGOUS
    is_masked = condition is ClinvarCondition.MASKED
    spiked = []
    for spike in case.variants:
        variants = annotated.get(spike.clinvar_id, [])
        if is_masked:
            variants = list(filter(_is_loaded_without_clinvar, variants))
        if not variants:
            return None
        for variant in variants:
            called = replace(
                variant,
                genotype="1/1" if is_homozygous else "0/1",
                filter_status="PASS",
                genotype_quality=99.0,
                depth=30.0,
                allele_depths=(0, 30) if is_homozygous else (15, 15),
                phase_set=None,
                phase_haplotype=None,
                original_site=None,
            )
            if is_masked:
                called = replace(
                    called,
                    clinvar_significance="",
                    clinvar_review_status="",
                    clinvar_conflict="",
                )
            spiked.append(called)
    return spiked


def _outcome(
    candidates: Sequence[Candidate], truth: frozenset[Allele], gene: str
) -> Outcome:
    ranked = list(enumerate(candidates, start=1))
    score = score_predictions(
        [
            (candidate.hypothesis.alleles, candidate.epcr)
            for candidate in candidates[: submission.MAX_ROWS]
        ],
        truth,
    )
    return Outcome(
        rank=next((rank for rank, c in ranked if c.hypothesis.alleles == truth), None),
        gene_rank=next((rank for rank, c in ranked if c.hypothesis.gene == gene), None),
        rank_points=score.rank_points,
        f_max=round(score.f_max, 4),
    )


def _phenotype_matches(
    model: PhenotypeModel, terms: Sequence[str] | None, genes: Collection[str]
) -> dict[str, PhenotypeMatch] | None:
    if terms is None:
        return None
    return {gene: model.score_gene(gene, terms) for gene in genes}


def evaluate(
    cases: Sequence[Case],
    annotated: Mapping[str, Sequence[Variant]],
    background: Sequence[Variant],
    gene_validity: Mapping[str, GeneValidity],
    model: PhenotypeModel,
    config: ScoringConfig,
) -> list[CaseResult]:
    """Rank spike-in cases under every ClinVar and phenotype condition.

    Skip cases whose causal alleles already occur in the background genome.
    """
    is_background_male = is_male(background)
    background_alleles = {variant.allele for variant in background}
    background_evidence = [
        Evidence.from_variant(variant, config) for variant in background
    ]
    gene_evidence: dict[str, list[Evidence]] = defaultdict(list)
    for item in background_evidence:
        gene_evidence[item.variant.gene].append(item)
    background_hypotheses = build_hypotheses(
        background_evidence, gene_validity, config, is_male=is_background_male
    )
    background_genes = {hypothesis.gene for hypothesis in background_hypotheses}

    results = []
    for case in cases:
        if any(variant.allele in background_alleles for variant in case.variants):
            continue
        visible = spiked_variants(case, annotated, ClinvarCondition.VISIBLE)
        if visible is None:
            results.append(
                CaseResult(
                    case,
                    dict.fromkeys(
                        itertools.product(ClinvarCondition, PhenotypeCondition),
                        NOT_RANKED,
                    ),
                )
            )
            continue
        if any(variant.allele in background_alleles for variant in visible):
            continue
        genes_to_score = background_genes | {variant.gene for variant in visible}
        phenotype_matches = {
            condition: _phenotype_matches(
                model, _terms_of_case(case, condition), genes_to_score
            )
            for condition in PhenotypeCondition
        }
        outcomes = {}
        for clinvar_condition in ClinvarCondition:
            spiked = spiked_variants(case, annotated, clinvar_condition)
            if spiked is None:
                for phenotype_condition in PhenotypeCondition:
                    outcomes[clinvar_condition, phenotype_condition] = NOT_RANKED
                continue
            truth = frozenset(variant.allele for variant in spiked)
            spiked_genes = sorted({variant.gene for variant in spiked})
            # Reuse background hypotheses for genes unaffected by the spike.
            evidence = [item for gene in spiked_genes for item in gene_evidence[gene]]
            evidence += [Evidence.from_variant(variant, config) for variant in spiked]
            hypotheses = [
                hypothesis
                for hypothesis in background_hypotheses
                if hypothesis.gene not in spiked_genes
            ]
            hypotheses += build_hypotheses(
                evidence, gene_validity, config, is_male=is_background_male
            )
            for phenotype_condition in PhenotypeCondition:
                candidates = rank_hypotheses(
                    hypotheses,
                    phenotype_matches=phenotype_matches[phenotype_condition],
                    config=config,
                )
                outcomes[clinvar_condition, phenotype_condition] = _outcome(
                    candidates, truth, case.gene
                )
        results.append(CaseResult(case, outcomes))
    return results


def write_case_table(results: Sequence[CaseResult], table_path: Path) -> None:
    conditions = list(itertools.product(ClinvarCondition, PhenotypeCondition))
    header = [
        "case",
        "type",
        "gene",
        "clinvar_ids",
        "kinds",
        "disease",
        "disease_matches_clinvar",
    ]
    for clinvar_condition, phenotype_condition in conditions:
        prefix = f"{phenotype_condition.value}_{clinvar_condition.value}"
        header += [
            f"{prefix}_{column}" for column in ("rank", "gene_rank", "points", "fmax")
        ]
    with table_path.open("w", encoding="utf-8", newline="") as table:
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for result in results:
            case = result.case
            row: list[object] = [
                case.name,
                case.genotype.value,
                case.gene,
                ";".join(variant.clinvar_id for variant in case.variants),
                case.kinds,
                case.disease,
                case.is_disease_named_by_clinvar,
            ]
            for condition in conditions:
                outcome = result.outcomes[condition]
                row += [
                    outcome.rank,
                    outcome.gene_rank,
                    outcome.rank_points,
                    outcome.f_max,
                ]
            writer.writerow(row)


def _percentage_ranked_within(outcomes: Sequence[Outcome], rank: int) -> float:
    ranked = sum(
        outcome.rank is not None and outcome.rank <= rank for outcome in outcomes
    )
    return round(100.0 * ranked / len(outcomes), 1)


def summarise(results: Sequence[CaseResult]) -> list[dict[str, str | int | float]]:
    summary: list[dict[str, str | int | float]] = []
    for clinvar_condition, phenotype_condition in itertools.product(
        ClinvarCondition, PhenotypeCondition
    ):
        for group, is_member in CASE_GROUPS:
            outcomes = [
                result.outcomes[clinvar_condition, phenotype_condition]
                for result in results
                if is_member(result.case)
            ]
            if not outcomes:
                continue
            count = len(outcomes)
            rank_points = sum(outcome.rank_points for outcome in outcomes)
            f_max = sum(outcome.f_max for outcome in outcomes)
            summary.append(
                {
                    "clinvar": clinvar_condition.value,
                    "phenotype": phenotype_condition.value,
                    "cases": group,
                    "n": count,
                    "top1_pct": _percentage_ranked_within(outcomes, 1),
                    "top3_pct": _percentage_ranked_within(outcomes, 3),
                    "top5_pct": _percentage_ranked_within(outcomes, 5),
                    "top10_pct": _percentage_ranked_within(outcomes, 10),
                    "mean_rank_points": round(rank_points / count, 1),
                    "mean_f_max": round(f_max / count, 3),
                }
            )
    return summary


def markdown_of_summary(summary: Sequence[Mapping[str, str | int | float]]) -> str:
    lines = [
        "| " + " | ".join(SUMMARY_COLUMNS) + " |",
        "|" + "---|" * len(SUMMARY_COLUMNS),
    ]
    lines += [
        "| " + " | ".join(str(row[column]) for column in SUMMARY_COLUMNS) + " |"
        for row in summary
    ]
    return "\n".join(lines) + "\n"
