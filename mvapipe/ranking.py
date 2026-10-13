import itertools
import math
import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from typing import NamedTuple, assert_never

from mvapipe.gene_validity import UNKNOWN_GENE, Gene_validity, Inheritance, Validity
from mvapipe.ontology import NO_MATCH, Phenotype_match
from mvapipe.scoring_config import (
    Call_quality,
    Frequency_window,
    Gene_prior,
    Impact_weights,
    Scoring_config,
)
from mvapipe.variant import Allele, Variant, Zygosity, alleles_of_genotype

LOSS_OF_FUNCTION = frozenset(
    {
        "stop_gained",
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "transcript_ablation",
    }
)
TRUNCATING = frozenset({"stop_gained", "frameshift_variant"})
INFRAME = frozenset({"inframe_insertion", "inframe_deletion", "protein_altering_variant"})
SPLICE_REGION = frozenset(
    {
        "splice_region_variant",
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
    }
)
# Substrings of ClinVar's review status, most reviewed first.
REVIEW_STARS = (
    ("practice_guideline", 4),
    ("reviewed_by_expert_panel", 3),
    ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
    ("criteria_provided,_conflicting", 1),
    ("criteria_provided,_single_submitter", 1),
)
MITOCHONDRIAL = ("MT", "M")
# GRCh38 positions between the two pseudo-autosomal regions.
X_NON_PSEUDO_AUTOSOMAL = range(2_781_480, 155_701_383)
MINIMUM_X_VARIANTS_TO_INFER_SEX = 20
MAXIMUM_X_HETEROZYGOSITY_OF_MALE = 0.15
RANKS_WITH_EPCR = 50
EPCR_CEILING = 0.99
EPCR_STEP = 1e-4


class Clinvar_class(Enum):
    Pathogenic = "pathogenic"
    Conflicting_with_pathogenic = "conflicting, some pathogenic"
    Conflicting_mostly_benign = "conflicting, mostly benign"
    Conflicting_some_benign = "conflicting, some benign"
    Conflicting = "conflicting"
    Benign = "benign"
    Uncertain = "uncertain"
    Unclassified = "unclassified"


class Model(Enum):
    Recessive_homozygous = "AR_hom"
    Recessive_compound = "AR_comphet"
    Recessive_single_hit = "AR_single"
    Dominant = "AD"
    X_linked = "XL"


class Phase(Enum):
    Cis = "cis"
    Trans = "trans"
    Unknown = "unknown"
    Not_applicable = "na"


class Conflict_counts(NamedTuple):
    pathogenic: int
    benign: int
    uncertain: int


@dataclass(frozen=True)
class Evidence:
    variant: Variant
    zygosity: Zygosity
    clinvar: Clinvar_class
    review_stars: int
    impact: float
    quality: float

    @staticmethod
    def of_variant(variant: Variant, config: Scoring_config) -> "Evidence":
        zygosity = Zygosity.of_genotype(variant.genotype)
        clinvar = clinvar_class(variant.clinvar_significance, variant.clinvar_conflict)
        stars = review_stars(variant.clinvar_review_status)
        return Evidence(
            variant=variant,
            zygosity=zygosity,
            clinvar=clinvar,
            review_stars=stars,
            impact=_impact(variant, clinvar, stars, config),
            quality=_quality(variant, zygosity, config.quality),
        )


@dataclass(frozen=True)
class Hypothesis:
    """One or two variants proposed as the cause under an inheritance model."""

    gene: str
    model: Model
    variants: tuple[Evidence] | tuple[Evidence, Evidence]
    genotype_evidence: float
    gene_prior: float
    phase: Phase = Phase.Not_applicable

    @property
    def alleles(self) -> frozenset[Allele]:
        return frozenset(evidence.variant.allele for evidence in self.variants)


@dataclass(frozen=True)
class Candidate:
    hypothesis: Hypothesis
    # None in a blind ranking.
    phenotype: Phenotype_match | None
    score: float
    # Decreases strictly with rank.
    epcr: float


def review_stars(review_status: str) -> int:
    review_status = review_status.lower()
    return next((stars for text, stars in REVIEW_STARS if text in review_status), 0)


def conflict_counts(conflict: str) -> Conflict_counts:
    """Count ClinVar submissions in fields such as 'Benign_(8)'."""
    counts = dict.fromkeys(Conflict_counts._fields, 0)
    for name, count in re.findall(r"([A-Za-z_/]+?)_?\((\d+)\)", conflict):
        for kind in counts:
            if kind in name.lower():
                counts[kind] += int(count)
    return Conflict_counts(**counts)


def clinvar_class(significance: str, conflict: str) -> Clinvar_class:
    significance = significance.lower()
    if "conflicting" in significance:
        counts = conflict_counts(conflict)
        if counts.pathogenic or "pathogenic" in conflict.lower():
            return Clinvar_class.Conflicting_with_pathogenic
        if counts.benign > counts.uncertain:
            return Clinvar_class.Conflicting_mostly_benign
        if counts.benign:
            return Clinvar_class.Conflicting_some_benign
        return Clinvar_class.Conflicting
    if "pathogenic" in significance:
        return Clinvar_class.Pathogenic
    if "benign" in significance:
        return Clinvar_class.Benign
    if "uncertain" in significance:
        return Clinvar_class.Uncertain
    return Clinvar_class.Unclassified


def _consequence_impact(variant: Variant, weights: Impact_weights) -> float:
    terms = set(variant.consequence.split("&"))
    if terms & LOSS_OF_FUNCTION:
        if terms & TRUNCATING and _is_in_last_exon(variant.exon):
            return weights.loss_of_function * weights.last_exon_truncation_factor
        return weights.loss_of_function
    if "start_lost" in terms:
        return weights.start_lost
    if "stop_lost" in terms:
        return weights.stop_lost
    if terms & INFRAME:
        return weights.inframe
    if "missense_variant" in terms:
        if variant.alphamissense is None:
            return weights.missense_without_score
        return weights.missense_base + weights.missense_slope * variant.alphamissense
    if terms & SPLICE_REGION:
        return weights.splice_region
    if "synonymous_variant" in terms:
        return weights.synonymous
    return weights.other


def _is_in_last_exon(exon: str) -> bool:
    """`exon` is VEP's 'number/total'; the number is a range for multi-exon deletions."""
    if "/" not in exon:
        return False
    number, total = exon.split("/")[:2]
    return number.split("-")[-1] == total


def _impact(
    variant: Variant, clinvar: Clinvar_class, stars: int, config: Scoring_config
) -> float:
    impact = _consequence_impact(variant, config.impact)
    weights = config.clinvar
    match clinvar:
        case Clinvar_class.Pathogenic:
            if stars >= 2:
                return max(impact, weights.pathogenic_two_stars_floor)
            if stars == 1:
                return max(impact, weights.pathogenic_one_star_floor)
            return max(impact, weights.pathogenic_no_stars_floor)
        case Clinvar_class.Conflicting_with_pathogenic:
            return max(impact, weights.conflicting_with_pathogenic_floor)
        case Clinvar_class.Benign:
            return min(impact, weights.benign_cap) if stars >= 1 else impact
        case Clinvar_class.Conflicting_mostly_benign:
            return min(impact, weights.conflicting_mostly_benign_cap)
        case Clinvar_class.Conflicting_some_benign:
            return impact * weights.conflicting_some_benign_factor
        case (
            Clinvar_class.Conflicting
            | Clinvar_class.Uncertain
            | Clinvar_class.Unclassified
        ):
            return impact
        case _ as unreachable:
            assert_never(unreachable)


def _quality(variant: Variant, zygosity: Zygosity, weights: Call_quality) -> float:
    quality = 1.0
    if variant.filter_status not in ("PASS", ".", ""):
        quality *= weights.failed_filter_factor
    genotype_quality = variant.genotype_quality
    if (
        genotype_quality is not None
        and genotype_quality < weights.minimum_genotype_quality
    ):
        quality *= weights.low_genotype_quality_factor
    allele_depths = variant.allele_depths
    depth = variant.depth
    if depth is None and allele_depths:
        depth = sum(allele_depths)
    if depth is not None and depth < weights.minimum_depth:
        quality *= weights.low_depth_factor
    is_heterozygous = zygosity is Zygosity.Heterozygous
    if is_heterozygous and allele_depths is not None and len(allele_depths) >= 2:
        total = sum(allele_depths)
        balance = allele_depths[-1] / total if total > 0 else None
        is_balanced = balance is None or (
            weights.minimum_allele_balance <= balance <= weights.maximum_allele_balance
        )
        if not is_balanced:
            quality *= weights.allele_balance_factor
    return quality


def _rarity_weight(allele_frequency: float, window: Frequency_window) -> float:
    if allele_frequency <= window.full:
        return 1.0
    if allele_frequency >= window.zero:
        return 0.0
    return (math.log10(window.zero) - math.log10(allele_frequency)) / (
        math.log10(window.zero) - math.log10(window.full)
    )


def _weight(evidence: Evidence, *, is_recessive: bool, config: Scoring_config) -> float:
    rarity = config.rarity
    allele_frequency = evidence.variant.allele_frequency
    is_exempt = (
        is_recessive
        and evidence.clinvar is Clinvar_class.Pathogenic
        and evidence.review_stars >= 2
        and allele_frequency < rarity.recessive_exempt_below
    )
    window = rarity.recessive if is_recessive else rarity.dominant
    rarity_weight = 1.0 if is_exempt else _rarity_weight(allele_frequency, window)
    return evidence.impact * rarity_weight * evidence.quality


def phase_of_pair(first: Variant, second: Variant) -> Phase:
    # Two alternate alleles of one site are on different haplotypes.
    if first.original_site and first.original_site == second.original_site:
        return Phase.Trans
    is_phased_together = (
        first.phase_set is not None
        and first.phase_set == second.phase_set
        and first.phase_haplotype is not None
        and second.phase_haplotype is not None
    )
    if not is_phased_together:
        return Phase.Unknown
    if first.phase_haplotype == second.phase_haplotype:
        return Phase.Cis
    return Phase.Trans


def _inheritance_of_model(model: Model) -> Inheritance:
    match model:
        case (
            Model.Recessive_homozygous
            | Model.Recessive_compound
            | Model.Recessive_single_hit
        ):
            return Inheritance.Autosomal_recessive
        case Model.Dominant:
            return Inheritance.Autosomal_dominant
        case Model.X_linked:
            return Inheritance.X_linked
        case _ as unreachable:
            assert_never(unreachable)


def _prior_without_evidence(
    inheritance: Inheritance, weights: Gene_prior, *, is_gene_unknown: bool
) -> float:
    match inheritance:
        case Inheritance.Autosomal_recessive:
            if is_gene_unknown:
                return weights.recessive_in_unknown_gene
            return weights.recessive_in_non_recessive_gene
        case Inheritance.Autosomal_dominant:
            if is_gene_unknown:
                return weights.dominant_in_unknown_gene
            return weights.dominant_in_non_dominant_gene
        case Inheritance.X_linked:
            return weights.x_linked_in_other_gene
        case _ as unreachable:
            assert_never(unreachable)


def _gene_prior(validity: Gene_validity, model: Model, weights: Gene_prior) -> float:
    inheritance = _inheritance_of_model(model)
    match validity.for_inheritance(inheritance):
        case Validity.Definitive:
            return weights.definitive
        case Validity.Strong:
            return weights.strong
        case Validity.Moderate:
            return weights.moderate
        case Validity.Supportive:
            return weights.supportive
        case Validity.Limited:
            return weights.limited
        case Validity.No_evidence:
            return _prior_without_evidence(
                inheritance,
                weights,
                is_gene_unknown=validity.best is Validity.No_evidence,
            )
        case _ as unreachable:
            assert_never(unreachable)


def _homozygous_hypotheses(
    gene: str,
    evidence: Sequence[Evidence],
    validity: Gene_validity,
    config: Scoring_config,
    *,
    is_on_x: bool,
    is_male: bool | None,
) -> list[Hypothesis]:
    hypotheses = []
    for item in evidence:
        if item.zygosity not in (Zygosity.Homozygous, Zygosity.Hemizygous):
            continue
        weight = _weight(item, is_recessive=True, config=config)
        if weight <= 0:
            continue
        is_x_linked = is_on_x and (is_male or item.zygosity is Zygosity.Hemizygous)
        model = Model.X_linked if is_x_linked else Model.Recessive_homozygous
        prior = _gene_prior(validity, model, config.gene_prior)
        hypotheses.append(Hypothesis(gene, model, (item,), weight, prior))
    return hypotheses


def _compound_hypotheses(
    gene: str,
    heterozygotes: Sequence[Evidence],
    validity: Gene_validity,
    config: Scoring_config,
) -> list[Hypothesis]:
    weights = config.inheritance
    weighted = sorted(
        (
            (_weight(item, is_recessive=True, config=config), item)
            for item in heterozygotes
        ),
        key=lambda pair: -pair[0],
    )
    strongest = [pair for pair in weighted if pair[0] >= weights.minimum_pair_weight]
    strongest = strongest[: weights.maximum_heterozygotes_per_gene]
    prior = _gene_prior(validity, Model.Recessive_compound, config.gene_prior)
    hypotheses = []
    for (first_weight, first), (second_weight, second) in itertools.combinations(
        strongest, 2
    ):
        phase = phase_of_pair(first.variant, second.variant)
        if phase is Phase.Cis:
            continue
        phase_factor = 1.0 if phase is Phase.Trans else weights.unknown_phase_factor
        genotype_evidence = math.sqrt(first_weight * second_weight) * phase_factor
        hypotheses.append(
            Hypothesis(
                gene,
                Model.Recessive_compound,
                (first, second),
                genotype_evidence,
                prior,
                phase,
            )
        )
    return hypotheses


def _single_hit_hypotheses(
    gene: str,
    heterozygotes: Sequence[Evidence],
    validity: Gene_validity,
    config: Scoring_config,
) -> list[Hypothesis]:
    weights = config.inheritance
    if validity.for_inheritance(Inheritance.Autosomal_recessive) < Validity.Supportive:
        return []
    prior = _gene_prior(validity, Model.Recessive_single_hit, config.gene_prior)
    hypotheses = []
    for item in heterozygotes:
        if item.impact < weights.single_hit_minimum_impact:
            continue
        genotype_evidence = weights.single_hit_factor * _weight(
            item, is_recessive=True, config=config
        )
        if genotype_evidence > 0:
            hypotheses.append(
                Hypothesis(
                    gene, Model.Recessive_single_hit, (item,), genotype_evidence, prior
                )
            )
    return hypotheses


def _dominant_hypotheses(
    gene: str,
    heterozygotes: Sequence[Evidence],
    validity: Gene_validity,
    config: Scoring_config,
    *,
    is_on_x: bool,
) -> list[Hypothesis]:
    weights = config.inheritance
    model = Model.X_linked if is_on_x else Model.Dominant
    prior = _gene_prior(validity, model, config.gene_prior)
    hypotheses = []
    for item in heterozygotes:
        if item.impact < weights.dominant_minimum_impact:
            continue
        weight = _weight(item, is_recessive=False, config=config)
        genotype_evidence = weight * weights.dominant_singleton_factor
        if genotype_evidence > 0:
            hypotheses.append(Hypothesis(gene, model, (item,), genotype_evidence, prior))
    return hypotheses


def hypotheses_of_evidence(
    evidence: Sequence[Evidence],
    gene_validity: Mapping[str, Gene_validity],
    config: Scoring_config,
    *,
    is_male: bool | None,
) -> list[Hypothesis]:
    gene_evidence: dict[str, list[Evidence]] = defaultdict(list)
    for item in evidence:
        if item.variant.allele.chromosome not in MITOCHONDRIAL:
            gene_evidence[item.variant.gene].append(item)

    hypotheses: list[Hypothesis] = []
    for gene, items in gene_evidence.items():
        validity = gene_validity.get(gene, UNKNOWN_GENE)
        is_on_x = items[0].variant.allele.chromosome == "X"
        heterozygotes = [item for item in items if item.zygosity is Zygosity.Heterozygous]
        hypotheses += _homozygous_hypotheses(
            gene, items, validity, config, is_on_x=is_on_x, is_male=is_male
        )
        if not is_on_x:
            hypotheses += _compound_hypotheses(gene, heterozygotes, validity, config)
            hypotheses += _single_hit_hypotheses(gene, heterozygotes, validity, config)
        hypotheses += _dominant_hypotheses(
            gene, heterozygotes, validity, config, is_on_x=is_on_x
        )
    return hypotheses


def _epcr_of_scores(scores: Sequence[float]) -> list[float]:
    """Decrease EPCR strictly for the first 50 candidates; use zero thereafter."""
    top_score = max(scores[0], 1e-9)
    ceiling = min(EPCR_CEILING, top_score)
    values: list[float] = []
    for score in scores[:RANKS_WITH_EPCR]:
        value = ceiling * (score / top_score) ** 2 if values else ceiling
        if values and value > values[-1] - EPCR_STEP:
            previous = values[-1]
            value = previous - EPCR_STEP if previous > 2 * EPCR_STEP else previous / 2
        values.append(value)
    return values + [0.0] * (len(scores) - len(values))


def rank_hypotheses(
    hypotheses: Sequence[Hypothesis],
    phenotype_matches: Mapping[str, Phenotype_match] | None,
    config: Scoring_config,
) -> list[Candidate]:
    """Keep the highest scoring hypothesis for each set of alleles."""
    floor = config.phenotype_floor
    scored: list[tuple[float, Hypothesis, Phenotype_match | None]] = []
    for hypothesis in hypotheses:
        phenotype = None
        phenotype_factor = 1.0
        if phenotype_matches is not None:
            phenotype = phenotype_matches.get(hypothesis.gene, NO_MATCH)
            phenotype_factor = floor + (1.0 - floor) * phenotype.similarity
        score = hypothesis.genotype_evidence * hypothesis.gene_prior * phenotype_factor
        scored.append((score, hypothesis, phenotype))
    scored.sort(
        key=lambda entry: (
            -entry[0],
            entry[1].gene,
            sorted(entry[1].alleles),
            entry[1].model.value,
        )
    )
    best: dict[frozenset[Allele], tuple[float, Hypothesis, Phenotype_match | None]] = {}
    for entry in scored:
        best.setdefault(entry[1].alleles, entry)
    ranked = list(best.values())
    if not ranked:
        return []
    epcr_values = _epcr_of_scores([score for score, _, _ in ranked])
    return [
        Candidate(hypothesis=hypothesis, phenotype=phenotype, score=score, epcr=epcr)
        for (score, hypothesis, phenotype), epcr in zip(ranked, epcr_values, strict=True)
    ]


def rank(
    variants: Sequence[Variant],
    gene_validity: Mapping[str, Gene_validity],
    phenotype_matches: Mapping[str, Phenotype_match] | None,
    config: Scoring_config,
) -> list[Candidate]:
    evidence = [Evidence.of_variant(variant, config) for variant in variants]
    hypotheses = hypotheses_of_evidence(
        evidence, gene_validity, config, is_male=is_male(variants)
    )
    return rank_hypotheses(hypotheses, phenotype_matches, config)


def is_male(variants: Sequence[Variant]) -> bool | None:
    """Infer sex from X heterozygosity; return None below the minimum variant count."""
    genotypes = [
        variant.genotype
        for variant in variants
        if variant.allele.chromosome == "X"
        and variant.allele.position in X_NON_PSEUDO_AUTOSOMAL
    ]
    if len(genotypes) < MINIMUM_X_VARIANTS_TO_INFER_SEX:
        return None
    heterozygous = sum(
        len(set(alleles_of_genotype(genotype))) > 1 for genotype in genotypes
    )
    return heterozygous / len(genotypes) < MAXIMUM_X_HETEROZYGOSITY_OF_MALE
