from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path

from benchmark.cases import AUTOSOMES, Case, SpikeVariant
from mvapipe.references import References
from mvapipe.stage_runner import StageRunner
from mvapipe.variant import Variant
from mvapipe.variant_table import load_variants
from mvapipe.vcf_steps import build_variant_table


def _genomic_order(variant: SpikeVariant) -> tuple[int, int, str, str]:
    chromosome, position, reference, alternate = variant.allele
    return AUTOSOMES.index(chromosome), position, reference, alternate


def write_spike_vcf(variants: Sequence[SpikeVariant], vcf_path: Path) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        *(f"##contig=<ID={chromosome}>" for chromosome in AUTOSOMES),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPIKE",
    ]
    for variant in sorted(variants, key=_genomic_order):
        chromosome, position, reference, alternate = variant.allele
        lines.append(
            f"{chromosome}\t{position}\t{variant.clinvar_id}\t{reference}\t{alternate}"
            "\t.\tPASS\t.\tGT\t0/1"
        )
    vcf_path.write_text("\n".join(lines) + "\n", encoding="utf8")


def annotate_spike_variants(
    cases: Sequence[Case], directory: Path, references: References, threads: int
) -> dict[str, list[Variant]]:
    """Annotate case variants and group them by ClinVar ID."""
    spike_variants = {
        variant.clinvar_id: variant for case in cases for variant in case.variants
    }
    directory.mkdir(parents=True, exist_ok=True)
    vcf_path = directory / "spikes.vcf"
    write_spike_vcf(list(spike_variants.values()), vcf_path)
    table_path = build_variant_table(
        vcf_path,
        sample=None,
        directory=directory,
        references=references,
        thread_count=threads,
        runner=StageRunner(directory),
    )
    variants = load_variants(table_path)
    variants_by_id: dict[str, list[Variant]] = defaultdict(list)
    for variant in variants:
        variants_by_id[variant.identifier].append(variant)
    return dict(variants_by_id)
