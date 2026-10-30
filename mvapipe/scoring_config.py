from dataclasses import dataclass


@dataclass(frozen=True)
class ImpactWeights:
    loss_of_function: float = 0.95
    start_lost: float = 0.70
    stop_lost: float = 0.60
    inframe: float = 0.50
    # Missense impact = base + slope * AlphaMissense score.
    missense_base: float = 0.15
    missense_slope: float = 0.65
    missense_without_score: float = 0.35
    splice_region: float = 0.25
    synonymous: float = 0.03
    other: float = 0.01
    last_exon_truncation_factor: float = 0.80


@dataclass(frozen=True)
class ClinvarWeights:
    pathogenic_two_stars_floor: float = 1.00
    pathogenic_one_star_floor: float = 0.95
    pathogenic_no_stars_floor: float = 0.85
    conflicting_with_pathogenic_floor: float = 0.60
    benign_cap: float = 0.02
    conflicting_mostly_benign_cap: float = 0.05
    conflicting_some_benign_factor: float = 0.50


@dataclass(frozen=True)
class FrequencyWindow:
    """Frequency thresholds for reducing the rarity weight from 1 to 0.

    The weight falls linearly in log frequency between these thresholds.
    """

    full_weight_frequency: float
    zero_weight_frequency: float


@dataclass(frozen=True)
class RarityWeights:
    recessive: FrequencyWindow = FrequencyWindow(
        full_weight_frequency=0.001, zero_weight_frequency=0.02
    )
    dominant: FrequencyWindow = FrequencyWindow(
        full_weight_frequency=0.00001, zero_weight_frequency=0.001
    )
    # Applies to ClinVar pathogenic alleles with two or more review stars.
    recessive_exempt_below: float = 0.05


@dataclass(frozen=True)
class CallQualityWeights:
    failed_filter_factor: float = 0.5
    minimum_genotype_quality: float = 20
    low_genotype_quality_factor: float = 0.5
    minimum_depth: float = 8
    low_depth_factor: float = 0.5
    minimum_allele_balance: float = 0.25
    maximum_allele_balance: float = 0.85
    allele_balance_factor: float = 0.7


@dataclass(frozen=True)
class InheritanceWeights:
    unknown_phase_factor: float = 0.90
    dominant_singleton_factor: float = 0.80
    dominant_minimum_impact: float = 0.15
    single_hit_factor: float = 0.25
    single_hit_minimum_impact: float = 0.85
    minimum_pair_weight: float = 0.05
    maximum_heterozygotes_per_gene: int = 12


@dataclass(frozen=True)
class GenePriorWeights:
    definitive: float = 1.0
    strong: float = 1.0
    moderate: float = 0.9
    supportive: float = 0.8
    limited: float = 0.55
    recessive_in_non_recessive_gene: float = 0.50
    recessive_in_unknown_gene: float = 0.45
    dominant_in_non_dominant_gene: float = 0.15
    dominant_in_unknown_gene: float = 0.30
    x_linked_in_other_gene: float = 0.40


@dataclass(frozen=True)
class ScoringConfig:
    impact: ImpactWeights = ImpactWeights()
    clinvar: ClinvarWeights = ClinvarWeights()
    rarity: RarityWeights = RarityWeights()
    quality: CallQualityWeights = CallQualityWeights()
    inheritance: InheritanceWeights = InheritanceWeights()
    gene_prior: GenePriorWeights = GenePriorWeights()
    phenotype_floor: float = 0.30
