from dataclasses import dataclass


@dataclass(frozen=True)
class Impact_weights:
    loss_of_function: float = 0.95
    start_lost: float = 0.70
    stop_lost: float = 0.60
    inframe: float = 0.50
    # missense = base + slope * AlphaMissense score
    missense_base: float = 0.15
    missense_slope: float = 0.65
    missense_without_score: float = 0.35
    splice_region: float = 0.25
    synonymous: float = 0.03
    other: float = 0.01
    last_exon_truncation_factor: float = 0.80


@dataclass(frozen=True)
class Clinvar_weights:
    pathogenic_two_stars_floor: float = 1.00
    pathogenic_one_star_floor: float = 0.95
    pathogenic_no_stars_floor: float = 0.85
    conflicting_with_pathogenic_floor: float = 0.60
    benign_cap: float = 0.02
    conflicting_mostly_benign_cap: float = 0.05
    conflicting_some_benign_factor: float = 0.50


@dataclass(frozen=True)
class Frequency_window:
    """Rarity weight falls from 1 to 0 between the bounds on a log scale."""

    full: float
    zero: float


@dataclass(frozen=True)
class Rarity:
    recessive: Frequency_window = Frequency_window(full=0.001, zero=0.02)
    dominant: Frequency_window = Frequency_window(full=0.00001, zero=0.001)
    # Applies to ClinVar pathogenic alleles with two or more review stars.
    recessive_exempt_below: float = 0.05


@dataclass(frozen=True)
class Call_quality:
    failed_filter_factor: float = 0.5
    minimum_genotype_quality: float = 20
    low_genotype_quality_factor: float = 0.5
    minimum_depth: float = 8
    low_depth_factor: float = 0.5
    minimum_allele_balance: float = 0.25
    maximum_allele_balance: float = 0.85
    allele_balance_factor: float = 0.7


@dataclass(frozen=True)
class Inheritance_weights:
    unknown_phase_factor: float = 0.90
    dominant_singleton_factor: float = 0.80
    dominant_minimum_impact: float = 0.15
    single_hit_factor: float = 0.25
    single_hit_minimum_impact: float = 0.85
    minimum_pair_weight: float = 0.05
    maximum_heterozygotes_per_gene: int = 12


@dataclass(frozen=True)
class Gene_prior:
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
class Scoring_config:
    impact: Impact_weights = Impact_weights()
    clinvar: Clinvar_weights = Clinvar_weights()
    rarity: Rarity = Rarity()
    quality: Call_quality = Call_quality()
    inheritance: Inheritance_weights = Inheritance_weights()
    gene_prior: Gene_prior = Gene_prior()
    phenotype_floor: float = 0.30
