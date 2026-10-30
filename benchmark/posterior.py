"""Evaluate complete autosomal genotypes with held-out genes and source exclusion."""

import csv
import random
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from enum import Enum
from pathlib import Path

from mvapipe.posterior.analysis import DEFAULT_CONFIG, PosteriorConfig, analyse_genome
from mvapipe.posterior.folds import assign_gene_folds
from mvapipe.posterior.hypotheses import RankedAlleleSets
from mvapipe.posterior.priors import Genotype
from mvapipe.posterior.training import SEED, TrainingCase, TrainingData
from mvapipe.variant import Allele, Variant


class PhenotypeCondition(Enum):
    REAL = "real"
    BLIND = "blind"
    NOISY = "noisy"


@dataclass(frozen=True)
class EvaluationOptions:
    """One evaluation condition; run conditions separately to keep calibration aligned."""

    phenotype: PhenotypeCondition = PhenotypeCondition.REAL
    mask_clinvar: bool = True
    remove_causes: bool = False
    posterior: PosteriorConfig = DEFAULT_CONFIG


@dataclass(frozen=True)
class PosteriorResult:
    case: str
    gene: str
    genotype: str
    fold: int
    status: str
    rank: int | None = None
    gene_rank: int | None = None
    truth_probability: float = 0.0
    unexplained: float | None = None
    top_probability: float | None = None


def spike_variants(
    case: TrainingCase,
    annotated: Mapping[str, Sequence[Variant]],
    *,
    mask_clinvar: bool,
) -> list[Variant] | None:
    """Create clean synthetic calls, returning None if a causal allele is unavailable.

    Masking removes ClinVar fields and annotations retained only because of ClinVar.
    Require an annotation for the reported causal gene for every allele.
    """
    is_homozygous = case.genotype is Genotype.HOMOZYGOUS
    spiked = []
    for key in case.allele_keys:
        variants = annotated.get(key, ())
        if mask_clinvar:
            variants = [
                variant
                for variant in variants
                if variant.vep_impact != "MODIFIER" and not variant.is_gene_from_clinvar
            ]
        if not any(variant.gene == case.gene for variant in variants):
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
            if mask_clinvar:
                called = replace(
                    called,
                    clinvar_significance="",
                    clinvar_review_status="",
                    clinvar_conflict="",
                )
            spiked.append(called)
    return spiked


def exclude_test_cases(training: TrainingData, case: TrainingCase) -> TrainingData:
    """Exclude the test gene's fold, case ID and overlapping normalised causal alleles."""
    folds = assign_gene_folds(
        [case.gene, *(item.gene for item in training.cases)],
        known_folds=training.gene_folds,
    )
    test_alleles = {
        variant.allele
        for key in case.allele_keys
        for variant in training.alleles.get(key, ())
    }
    kept = [
        item
        for item in training.cases
        if item.identifier != case.identifier
        and folds[item.gene] != folds[case.gene]
        and not test_alleles.intersection(
            variant.allele
            for key in item.allele_keys
            for variant in training.alleles.get(key, ())
        )
    ]
    return replace(training, cases=kept)


def _patient_terms(
    case: TrainingCase, training: TrainingData, condition: PhenotypeCondition
) -> tuple[str, ...] | None:
    if condition is PhenotypeCondition.BLIND:
        return None
    if condition is PhenotypeCondition.REAL:
        return case.terms
    generator = random.Random(f"{SEED}:{case.identifier}")
    ontology = training.similarity.ontology
    kept = [term for term in case.terms if generator.random() >= 0.5]
    if not kept:
        kept = [generator.choice(case.terms)]
    noisy = []
    for term in kept:
        parents = sorted(
            p for p in ontology.parents.get(term, ()) if ontology.is_feature(p)
        )
        noisy.append(
            generator.choice(parents) if parents and generator.random() < 0.5 else term
        )
    annotation_terms = sorted(
        {
            term
            for diseases in training.similarity.disease_terms.values()
            for terms in diseases.values()
            for term in terms
        }
    )
    noisy.extend(generator.sample(annotation_terms, 2))
    return tuple(dict.fromkeys(noisy))


def outcome(
    ranked: RankedAlleleSets, truth: frozenset[Allele], gene: str
) -> tuple[int | None, int | None, float]:
    """Return conservative genotype and gene ranks and the true label's probability.

    The genotype rank includes the unexplained label. The gene rank sums all model
    mass for each gene and excludes the unexplained label. Ties receive the worst rank.
    An empty truth denotes a cause-removal control whose correct label is unexplained.
    """
    probabilities = dict(zip(ranked.allele_sets, ranked.posteriors.tolist(), strict=True))
    probabilities[frozenset()] = ranked.unexplained
    probability = probabilities.get(truth, 0.0)
    rank = (
        sum(value >= probability for value in probabilities.values())
        if truth in probabilities
        else None
    )
    gene_probability = ranked.gene_posteriors.get(gene)
    gene_rank = (
        sum(value >= gene_probability for value in ranked.gene_posteriors.values())
        if gene_probability is not None
        else None
    )
    return rank, gene_rank, probability


def evaluate_cases(
    cases: Sequence[TrainingCase],
    background: Sequence[Variant],
    training: TrainingData,
    *,
    options: EvaluationOptions,
) -> list[PosteriorResult]:
    """Refit on each augmented genome after excluding the test fold from training.

    Keep one result per supplied case, including unannotated and unsupported cases.
    Missing PubMed IDs prevent source-excluded phenotype evaluation. Fitting failures
    propagate. The background is a fixed genome; synthetic calls do not model patients.
    """
    identifiers = [case.identifier for case in cases]
    if len(set(identifiers)) != len(identifiers):
        message = "evaluation case IDs must be unique"
        raise ValueError(message)
    folds = assign_gene_folds(
        (case.gene for case in cases), known_folds=training.gene_folds
    )
    results = []
    for case in cases:
        result = PosteriorResult(
            case.identifier, case.gene, case.genotype.value, folds[case.gene], "ranked"
        )
        if case.genotype is Genotype.SINGLE_RECESSIVE:
            results.append(replace(result, status="incomplete_genotype"))
            continue
        if options.phenotype is not PhenotypeCondition.BLIND and not case.pubmed_ids:
            results.append(replace(result, status="missing_publications"))
            continue
        spiked = spike_variants(case, training.alleles, mask_clinvar=options.mask_clinvar)
        if spiked is None:
            results.append(replace(result, status="missing_annotation"))
            continue
        truth = frozenset(variant.allele for variant in spiked)
        # Replace existing calls at the causal alleles to avoid duplicate gene rows.
        genome = [variant for variant in background if variant.allele not in truth]
        if options.remove_causes:
            truth = frozenset()
        else:
            genome.extend(spiked)
        analysis = analyse_genome(
            genome,
            _patient_terms(case, training, options.phenotype),
            exclude_test_cases(training, case),
            top=None,
            config=options.posterior,
            patient_publications=case.pubmed_ids,
        )
        rank, gene_rank, probability = outcome(analysis.ranked, truth, case.gene)
        results.append(
            replace(
                result,
                rank=rank,
                gene_rank=gene_rank,
                truth_probability=probability,
                unexplained=analysis.ranked.unexplained,
                top_probability=max(
                    analysis.ranked.unexplained,
                    max(analysis.ranked.posteriors.tolist(), default=0.0),
                ),
            )
        )
    return results


def write_results(results: Sequence[PosteriorResult], path: Path) -> None:
    """Write all cases, retaining failed records in the evaluation denominator."""
    with path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(
            destination, fieldnames=PosteriorResult.__dataclass_fields__
        )
        writer.writeheader()
        writer.writerows(asdict(result) for result in results)
