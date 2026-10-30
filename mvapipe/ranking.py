import itertools
import math
import re
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from typing import NamedTuple, assert_never

from mvapipe.gene_validity import UNKNOWN_GENE, GeneValidity, Inheritance, Validity
from mvapipe.ontology import NO_MATCH, PhenotypeMatch
from mvapipe.scoring_config import (
    CallQualityWeights,
    FrequencyWindow,
    GenePriorWeights,
    ImpactWeights,
    ScoringConfig,
)
from mvapipe.variant import Allele, Variant, Zygosity, split_genotype

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
# Match ClinVar review-status text from highest to lowest star rating.
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
X_PSEUDO_AUTOSOMAL = (range(10_001, 2_781_480), range(155_701_383, 156_030_896))
MINIMUM_X_VARIANTS_TO_INFER_SEX = 20
MAXIMUM_X_HETEROZYGOSITY_OF_MALE = 0.15
RANKS_WITH_EPCR = 50
EPCR_CEILING = 0.99
EPCR_STEP = 1e-4


class ClinvarClass(Enum):
    PATHOGENIC = "pathogenic"
    CONFLICTING_WITH_PATHOGENIC = "conflicting, some pathogenic"
    CONFLICTING_MOSTLY_BENIGN = "conflicting, mostly benign"
    CONFLICTING_SOME_BENIGN = "conflicting, some benign"
    CONFLICTING = "conflicting"
    BENIGN = "benign"
    UNCERTAIN = "uncertain"
    UNCLASSIFIED = "unclassified"


class Model(Enum):
    RECESSIVE_HOMOZYGOUS = "AR_hom"
    RECESSIVE_COMPOUND = "AR_comphet"
    RECESSIVE_SINGLE_HIT = "AR_single"
    DOMINANT = "AD"
    X_LINKED = "XL"


class Phase(Enum):
    CIS = "cis"
    TRANS = "trans"
    UNKNOWN = "unknown"
    NOT_APPLICABLE = "na"


class ConflictCounts(NamedTuple):
    pathogenic: int
    benign: int
    uncertain: int


@dataclass(frozen=True)
class Evidence:
    variant: Variant
    zygosity: Zygosity
    clinvar: ClinvarClass
    review_stars: int
    impact: float
    quality: float

    @staticmethod
    def from_variant(variant: Variant, config: ScoringConfig) -> "Evidence":
        zygosity = Zygosity.from_genotype(variant.genotype)
        clinvar = classify_clinvar(
            significance=variant.clinvar_significance, conflict=variant.clinvar_conflict
        )
        stars = parse_review_stars(variant.clinvar_review_status)
        return Evidence(
            variant=variant,
            zygosity=zygosity,
            clinvar=clinvar,
            review_stars=stars,
            impact=_compute_impact(variant, clinvar, stars, config),
            quality=_compute_quality(variant, zygosity, config.quality),
        )


@dataclass(frozen=True, kw_only=True)
class Hypothesis:
    """One or two candidate causal variants under a single inheritance model."""

    gene: str
    model: Model
    evidence: tuple[Evidence] | tuple[Evidence, Evidence]
    genotype_evidence: float
    gene_prior: float
    phase: Phase = Phase.NOT_APPLICABLE

    @property
    def alleles(self) -> frozenset[Allele]:
        return frozenset(evidence.variant.allele for evidence in self.evidence)


@dataclass(frozen=True)
class Candidate:
    hypothesis: Hypothesis
    phenotype: PhenotypeMatch | None
    score: float
    epcr: float


def parse_review_stars(review_status: str) -> int:
    """Convert a ClinVar review status to a star rating from 0 to 4."""
    review_status = review_status.lower()
    return next((stars for text, stars in REVIEW_STARS if text in review_status), 0)


def parse_conflict_counts(conflict: str) -> ConflictCounts:
    """Count ClinVar submissions in fields such as 'Benign_(8)'."""
    counts = dict.fromkeys(ConflictCounts._fields, 0)
    for name, count in re.findall(r"([A-Za-z_/]+?)_?\((\d+)\)", conflict):
        for kind in counts:
            if kind in name.lower():
                counts[kind] += int(count)
    return ConflictCounts(**counts)


def classify_clinvar(*, significance: str, conflict: str) -> ClinvarClass:
    """Classify CLNSIG and summarise conflicting submissions from CLNSIGCONF."""
    significance = significance.lower()
    if "conflicting" in significance:
        counts = parse_conflict_counts(conflict)
        if counts.pathogenic or "pathogenic" in conflict.lower():
            return ClinvarClass.CONFLICTING_WITH_PATHOGENIC
        if counts.benign > counts.uncertain:
            return ClinvarClass.CONFLICTING_MOSTLY_BENIGN
        if counts.benign:
            return ClinvarClass.CONFLICTING_SOME_BENIGN
        return ClinvarClass.CONFLICTING
    if "pathogenic" in significance:
        return ClinvarClass.PATHOGENIC
    if "benign" in significance:
        return ClinvarClass.BENIGN
    if "uncertain" in significance:
        return ClinvarClass.UNCERTAIN
    return ClinvarClass.UNCLASSIFIED


def _compute_consequence_impact(variant: Variant, weights: ImpactWeights) -> float:
    terms = set(variant.consequence.split("&"))
    if terms & LOSS_OF_FUNCTION:
        if terms & TRUNCATING and is_in_last_exon(variant.exon):
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


def is_in_last_exon(exon: str) -> bool:
    """Check whether VEP's 'number/total' exon field reaches the last exon.

    For a range such as '3-5/5', compare the range end with the total.
    """
    if "/" not in exon:
        return False
    number, total = exon.split("/")[:2]
    return number.split("-")[-1] == total


def _compute_impact(
    variant: Variant, clinvar: ClinvarClass, stars: int, config: ScoringConfig
) -> float:
    impact = _compute_consequence_impact(variant, config.impact)
    weights = config.clinvar
    match clinvar:
        case ClinvarClass.PATHOGENIC:
            if stars >= 2:
                return max(impact, weights.pathogenic_two_stars_floor)
            if stars == 1:
                return max(impact, weights.pathogenic_one_star_floor)
            return max(impact, weights.pathogenic_no_stars_floor)
        case ClinvarClass.CONFLICTING_WITH_PATHOGENIC:
            return max(impact, weights.conflicting_with_pathogenic_floor)
        case ClinvarClass.BENIGN:
            return min(impact, weights.benign_cap) if stars >= 1 else impact
        case ClinvarClass.CONFLICTING_MOSTLY_BENIGN:
            return min(impact, weights.conflicting_mostly_benign_cap)
        case ClinvarClass.CONFLICTING_SOME_BENIGN:
            return impact * weights.conflicting_some_benign_factor
        case (
            ClinvarClass.CONFLICTING | ClinvarClass.UNCERTAIN | ClinvarClass.UNCLASSIFIED
        ):
            return impact
        case _ as unreachable:
            assert_never(unreachable)


def _compute_quality(
    variant: Variant, zygosity: Zygosity, weights: CallQualityWeights
) -> float:
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
    is_heterozygous = zygosity is Zygosity.HETEROZYGOUS
    if is_heterozygous and allele_depths is not None and len(allele_depths) >= 2:
        total = sum(allele_depths)
        balance = allele_depths[-1] / total if total > 0 else None
        is_balanced = balance is None or (
            weights.minimum_allele_balance <= balance <= weights.maximum_allele_balance
        )
        if not is_balanced:
            quality *= weights.allele_balance_factor
    return quality


def _compute_rarity_weight(
    allele_frequency: float | None, window: FrequencyWindow
) -> float:
    """Compute the rarity weight, using 1 when the allele frequency is missing."""
    if allele_frequency is None or allele_frequency <= window.full_weight_frequency:
        return 1.0
    if allele_frequency >= window.zero_weight_frequency:
        return 0.0
    return (math.log10(window.zero_weight_frequency) - math.log10(allele_frequency)) / (
        math.log10(window.zero_weight_frequency)
        - math.log10(window.full_weight_frequency)
    )


def _compute_weight(
    evidence: Evidence, *, is_recessive: bool, config: ScoringConfig
) -> float:
    rarity = config.rarity
    allele_frequency = evidence.variant.allele_frequency
    is_exempt = (
        is_recessive
        and evidence.clinvar is ClinvarClass.PATHOGENIC
        and evidence.review_stars >= 2
        and allele_frequency is not None
        and allele_frequency < rarity.recessive_exempt_below
    )
    window = rarity.recessive if is_recessive else rarity.dominant
    rarity_weight = 1.0 if is_exempt else _compute_rarity_weight(allele_frequency, window)
    return evidence.impact * rarity_weight * evidence.quality


def infer_phase(first: Variant, second: Variant) -> Phase:
    """Infer cis or trans phase for two heterozygous variants in the same sample.

    Return UNKNOWN unless a shared original site or phase set establishes the phase.
    """
    # Two alternate alleles of one site are on different haplotypes.
    if first.original_site and first.original_site == second.original_site:
        return Phase.TRANS
    is_phased_together = (
        first.phase_set is not None
        and first.phase_set == second.phase_set
        and first.phase_haplotype is not None
        and second.phase_haplotype is not None
    )
    if not is_phased_together:
        return Phase.UNKNOWN
    if first.phase_haplotype == second.phase_haplotype:
        return Phase.CIS
    return Phase.TRANS


def _map_model_to_inheritance(model: Model) -> Inheritance:
    match model:
        case (
            Model.RECESSIVE_HOMOZYGOUS
            | Model.RECESSIVE_COMPOUND
            | Model.RECESSIVE_SINGLE_HIT
        ):
            return Inheritance.AUTOSOMAL_RECESSIVE
        case Model.DOMINANT:
            return Inheritance.AUTOSOMAL_DOMINANT
        case Model.X_LINKED:
            return Inheritance.X_LINKED
        case _ as unreachable:
            assert_never(unreachable)


def _compute_prior_without_evidence(
    inheritance: Inheritance, weights: GenePriorWeights, *, is_gene_unknown: bool
) -> float:
    match inheritance:
        case Inheritance.AUTOSOMAL_RECESSIVE:
            if is_gene_unknown:
                return weights.recessive_in_unknown_gene
            return weights.recessive_in_non_recessive_gene
        case Inheritance.AUTOSOMAL_DOMINANT:
            if is_gene_unknown:
                return weights.dominant_in_unknown_gene
            return weights.dominant_in_non_dominant_gene
        case Inheritance.X_LINKED:
            return weights.x_linked_in_other_gene
        case _ as unreachable:
            assert_never(unreachable)


def _compute_gene_prior(
    validity: GeneValidity, model: Model, weights: GenePriorWeights
) -> float:
    inheritance = _map_model_to_inheritance(model)
    match validity.for_inheritance(inheritance):
        case Validity.DEFINITIVE:
            return weights.definitive
        case Validity.STRONG:
            return weights.strong
        case Validity.MODERATE:
            return weights.moderate
        case Validity.SUPPORTIVE:
            return weights.supportive
        case Validity.LIMITED:
            return weights.limited
        case Validity.NO_EVIDENCE:
            return _compute_prior_without_evidence(
                inheritance,
                weights,
                is_gene_unknown=validity.best is Validity.NO_EVIDENCE,
            )
        case _ as unreachable:
            assert_never(unreachable)


def _build_homozygous_hypotheses(
    gene: str,
    evidence: Sequence[Evidence],
    validity: GeneValidity,
    config: ScoringConfig,
    *,
    is_on_x: bool,
    is_male: bool | None,
) -> list[Hypothesis]:
    hypotheses = []
    for item in evidence:
        if item.zygosity not in (Zygosity.HOMOZYGOUS, Zygosity.HEMIZYGOUS):
            continue
        weight = _compute_weight(item, is_recessive=True, config=config)
        if weight <= 0:
            continue
        is_x_linked = is_on_x and (is_male or item.zygosity is Zygosity.HEMIZYGOUS)
        model = Model.X_LINKED if is_x_linked else Model.RECESSIVE_HOMOZYGOUS
        prior = _compute_gene_prior(validity, model, config.gene_prior)
        hypotheses.append(
            Hypothesis(
                gene=gene,
                model=model,
                evidence=(item,),
                genotype_evidence=weight,
                gene_prior=prior,
            )
        )
    return hypotheses


def _build_compound_hypotheses(
    gene: str,
    heterozygotes: Sequence[Evidence],
    validity: GeneValidity,
    config: ScoringConfig,
) -> list[Hypothesis]:
    weights = config.inheritance
    weighted = sorted(
        (
            (_compute_weight(item, is_recessive=True, config=config), item)
            for item in heterozygotes
        ),
        key=lambda pair: -pair[0],
    )
    strongest = [pair for pair in weighted if pair[0] >= weights.minimum_pair_weight]
    strongest = strongest[: weights.maximum_heterozygotes_per_gene]
    prior = _compute_gene_prior(validity, Model.RECESSIVE_COMPOUND, config.gene_prior)
    hypotheses = []
    for (first_weight, first), (second_weight, second) in itertools.combinations(
        strongest, 2
    ):
        if first.variant.allele == second.variant.allele:
            continue
        phase = infer_phase(first.variant, second.variant)
        if phase is Phase.CIS:
            continue
        phase_factor = 1.0 if phase is Phase.TRANS else weights.unknown_phase_factor
        genotype_evidence = math.sqrt(first_weight * second_weight) * phase_factor
        hypotheses.append(
            Hypothesis(
                gene=gene,
                model=Model.RECESSIVE_COMPOUND,
                evidence=(first, second),
                genotype_evidence=genotype_evidence,
                gene_prior=prior,
                phase=phase,
            )
        )
    return hypotheses


def _build_single_hit_hypotheses(
    gene: str,
    heterozygotes: Sequence[Evidence],
    validity: GeneValidity,
    config: ScoringConfig,
) -> list[Hypothesis]:
    weights = config.inheritance
    if validity.for_inheritance(Inheritance.AUTOSOMAL_RECESSIVE) < Validity.SUPPORTIVE:
        return []
    prior = _compute_gene_prior(validity, Model.RECESSIVE_SINGLE_HIT, config.gene_prior)
    hypotheses = []
    for item in heterozygotes:
        if item.impact < weights.single_hit_minimum_impact:
            continue
        genotype_evidence = weights.single_hit_factor * _compute_weight(
            item, is_recessive=True, config=config
        )
        if genotype_evidence > 0:
            hypotheses.append(
                Hypothesis(
                    gene=gene,
                    model=Model.RECESSIVE_SINGLE_HIT,
                    evidence=(item,),
                    genotype_evidence=genotype_evidence,
                    gene_prior=prior,
                )
            )
    return hypotheses


def _build_dominant_hypotheses(
    gene: str,
    heterozygotes: Sequence[Evidence],
    validity: GeneValidity,
    config: ScoringConfig,
    *,
    is_on_x: bool,
) -> list[Hypothesis]:
    weights = config.inheritance
    model = Model.X_LINKED if is_on_x else Model.DOMINANT
    prior = _compute_gene_prior(validity, model, config.gene_prior)
    hypotheses = []
    for item in heterozygotes:
        if item.impact < weights.dominant_minimum_impact:
            continue
        weight = _compute_weight(item, is_recessive=False, config=config)
        genotype_evidence = weight * weights.dominant_singleton_factor
        if genotype_evidence > 0:
            hypotheses.append(
                Hypothesis(
                    gene=gene,
                    model=model,
                    evidence=(item,),
                    genotype_evidence=genotype_evidence,
                    gene_prior=prior,
                )
            )
    return hypotheses


def build_hypotheses(
    evidence: Sequence[Evidence],
    gene_validity: Mapping[str, GeneValidity],
    config: ScoringConfig,
    *,
    is_male: bool | None,
) -> list[Hypothesis]:
    """Build each gene's inheritance hypotheses, excluding mitochondrial variants."""
    gene_evidence: dict[str, list[Evidence]] = defaultdict(list)
    for item in evidence:
        if item.variant.allele.chromosome not in MITOCHONDRIAL:
            gene_evidence[item.variant.gene].append(item)

    hypotheses: list[Hypothesis] = []
    for gene, items in gene_evidence.items():
        validity = gene_validity.get(gene, UNKNOWN_GENE)
        allele = items[0].variant.allele
        is_on_x = allele.chromosome == "X" and not any(
            allele.position in region for region in X_PSEUDO_AUTOSOMAL
        )
        heterozygotes = [item for item in items if item.zygosity is Zygosity.HETEROZYGOUS]
        hypotheses += _build_homozygous_hypotheses(
            gene, items, validity, config, is_on_x=is_on_x, is_male=is_male
        )
        if not is_on_x:
            hypotheses += _build_compound_hypotheses(
                gene, heterozygotes, validity, config
            )
            hypotheses += _build_single_hit_hypotheses(
                gene, heterozygotes, validity, config
            )
        hypotheses += _build_dominant_hypotheses(
            gene, heterozygotes, validity, config, is_on_x=is_on_x
        )
    return hypotheses


def _assign_epcr(scores: Sequence[float]) -> list[float]:
    """Assign EPCR values to scores supplied in descending order.

    Values decrease strictly for the first RANKS_WITH_EPCR scores and are zero thereafter.
    Expect a nonempty sequence; raise ValueError if the top score is not positive.
    """
    top_score = scores[0]
    if top_score <= 0:
        message = f"the top score is {top_score}; EPCR needs a positive score"
        raise ValueError(message)
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
    *,
    phenotype_matches: Mapping[str, PhenotypeMatch] | None,
    config: ScoringConfig,
) -> list[Candidate]:
    """Score hypotheses and keep the highest-scoring one for each allele set.

    Sort by descending score, breaking ties by gene, alleles and model. If
    `phenotype_matches` is None, use a phenotype factor of 1 and leave matches unset.
    Otherwise, genes missing from the mapping have similarity 0.
    Return an empty list for no hypotheses; raise ValueError for a nonpositive top score.
    """
    floor = config.phenotype_floor
    scored: list[tuple[float, Hypothesis, PhenotypeMatch | None]] = []
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
    best: dict[frozenset[Allele], tuple[float, Hypothesis, PhenotypeMatch | None]] = {}
    for entry in scored:
        best.setdefault(entry[1].alleles, entry)
    ranked = list(best.values())
    if not ranked:
        return []
    epcr_values = _assign_epcr([score for score, _, _ in ranked])
    return [
        Candidate(hypothesis=hypothesis, phenotype=phenotype, score=score, epcr=epcr)
        for (score, hypothesis, phenotype), epcr in zip(ranked, epcr_values, strict=True)
    ]


def rank(
    variants: Sequence[Variant],
    gene_validity: Mapping[str, GeneValidity],
    *,
    phenotype_matches: Mapping[str, PhenotypeMatch] | None,
    config: ScoringConfig,
) -> list[Candidate]:
    """Build and rank variant hypotheses using `rank_hypotheses`."""
    evidence = [Evidence.from_variant(variant, config) for variant in variants]
    hypotheses = build_hypotheses(
        evidence, gene_validity, config, is_male=is_male(variants)
    )
    return rank_hypotheses(hypotheses, phenotype_matches=phenotype_matches, config=config)


def is_male(variants: Sequence[Variant]) -> bool | None:
    """Infer male status from heterozygosity outside the X pseudo-autosomal regions.

    Return None when fewer than MINIMUM_X_VARIANTS_TO_INFER_SEX calls are available.
    """
    genotypes = [
        variant.genotype
        for variant in variants
        if variant.allele.chromosome == "X"
        and variant.allele.position in X_NON_PSEUDO_AUTOSOMAL
    ]
    if len(genotypes) < MINIMUM_X_VARIANTS_TO_INFER_SEX:
        return None
    heterozygous_count = sum(
        len(set(split_genotype(genotype))) > 1 for genotype in genotypes
    )
    return heterozygous_count / len(genotypes) < MAXIMUM_X_HETEROZYGOSITY_OF_MALE
