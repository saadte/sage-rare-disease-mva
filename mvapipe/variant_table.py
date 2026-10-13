import csv
import gzip
import itertools
import re
from collections import Counter
from pathlib import Path

from mvapipe.variant import Allele, Variant, Zygosity, alleles_of_genotype

SITE_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "FILTER"]
INFO_COLUMNS = ["CLNSIG", "CLNREVSTAT", "CLNSIGCONF", "GENEINFO", "ORIG"]
OPTIONAL_FORMAT_COLUMNS = ["AD", "DP", "GQ", "PGT", "PID", "PS"]
CONSEQUENCE_COLUMNS = [
    "SYMBOL",
    "Consequence",
    "IMPACT",
    "MANE_SELECT",
    "EXON",
    "Protein_position",
    "Amino_acids",
    "gnomADe_AF",
    "gnomADg_AF",
    "MAX_AF",
    "am_pathogenicity",
]

# VEP's IMPACT values, most severe first.
IMPACTS = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

Row = dict[str, str]


def load_variants(table_path: Path) -> tuple[list[Variant], Counter[str]]:
    """Read coding and ClinVar-flagged variants, one per allele and gene."""
    counts: Counter[str] = Counter()
    variants: list[Variant] = []
    with gzip.open(table_path, "rt", encoding="utf8", newline="") as table:
        rows = csv.DictReader(table, delimiter="\t")
        if rows.fieldnames is None:
            message = f"{table_path} has no header"
            raise ValueError(message)
        for _, rows_of_allele in itertools.groupby(rows, key=_allele_of_row):
            variants.extend(_variants_of_allele(list(rows_of_allele), counts))
    counts["kept"] = len(variants)
    return variants, counts


def _variants_of_allele(rows: list[Row], counts: Counter[str]) -> list[Variant]:
    counts["rows"] += len(rows)
    is_flagged = is_flagged_by_clinvar(rows[0])
    clinvar_genes = _clinvar_genes(rows[0]) if is_flagged else []
    # Assign ClinVar evidence to its named genes when VEP has rows for them.
    has_clinvar_gene = any(row["SYMBOL"] in clinvar_genes for row in rows)

    by_gene: dict[str, Variant] = {}
    # VEP writes a variant's rows in an order that changes from run to run.
    for row in sorted(rows, key=_severity_order):
        symbol = row["SYMBOL"]
        is_other_gene = (
            has_clinvar_gene and not _is_missing(symbol) and symbol not in clinvar_genes
        )
        if is_other_gene:
            counts["clinvar_other_gene"] += 1
        is_flagged_here = is_flagged and not is_other_gene

        if not _is_missing(symbol):
            gene, is_gene_from_clinvar = symbol, False
        elif is_flagged_here and clinvar_genes:
            gene, is_gene_from_clinvar = clinvar_genes[0], True
            counts["gene_from_clinvar"] += 1
        else:
            counts["no_gene"] += 1
            continue
        if row["IMPACT"] == "MODIFIER" and not is_flagged_here:
            counts["modifier_dropped"] += 1
            continue

        # Prefer VEP's gene assignment over a gene inferred from ClinVar.
        existing = by_gene.get(gene)
        if existing is not None:
            counts["duplicate_gene_rows"] += 1
            if is_gene_from_clinvar or not existing.is_gene_from_clinvar:
                continue
        by_gene[gene] = _variant_of_row(
            row,
            gene=gene,
            is_gene_from_clinvar=is_gene_from_clinvar,
            keep_clinvar=not is_other_gene,
        )
    return list(by_gene.values())


def _severity_order(row: Row) -> tuple[int, str, str, str, str]:
    return (
        IMPACTS.index(row["IMPACT"]),
        row["SYMBOL"],
        row["Consequence"],
        row["EXON"],
        row["am_pathogenicity"],
    )


def _variant_of_row(
    row: Row, *, gene: str, is_gene_from_clinvar: bool, keep_clinvar: bool
) -> Variant:
    frequencies = [
        frequency
        for column in ("gnomADe_AF", "gnomADg_AF", "MAX_AF")
        if (frequency := _largest_number(row[column])) is not None
    ]
    phase_set, phase_haplotype = _phase_of_row(row)
    return Variant(
        allele=_allele_of_row(row),
        gene=gene,
        consequence=row["Consequence"],
        vep_impact=row["IMPACT"],
        exon=_text(row["EXON"]),
        genotype=row["GT"],
        filter_status=row["FILTER"],
        genotype_quality=_largest_number(row.get("GQ", ".")),
        depth=_largest_number(row.get("DP", ".")),
        allele_depths=_allele_depths(row.get("AD", ".")),
        phase_set=phase_set,
        phase_haplotype=phase_haplotype,
        original_site=_original_site(row["ORIG"]),
        allele_frequency=max(frequencies, default=0.0),
        alphamissense=_largest_number(row["am_pathogenicity"]),
        clinvar_significance=_text(row["CLNSIG"]) if keep_clinvar else "",
        clinvar_review_status=_text(row["CLNREVSTAT"]) if keep_clinvar else "",
        clinvar_conflict=_text(row["CLNSIGCONF"]) if keep_clinvar else "",
        identifier=row["ID"],
        is_gene_from_clinvar=is_gene_from_clinvar,
    )


def _phase_of_row(row: Row) -> tuple[str | None, int | None]:
    """Return the phase set and the alternate allele's haplotype."""
    genotype = row["GT"]
    haplotype = _alternate_haplotype(genotype)
    physical_genotype, physical_set = row.get("PGT", "."), row.get("PID", ".")
    if not _is_missing(physical_set) and not _is_missing(physical_genotype):
        haplotypes = physical_genotype.split("|")
        if "1" in haplotypes:
            haplotype = haplotypes.index("1")
        return f"PID:{physical_set}", haplotype
    phase_set = row.get("PS", ".")
    if not _is_missing(phase_set) and "|" in genotype:
        return f"PS:{row['CHROM']}:{phase_set}", haplotype
    return None, haplotype


def _alternate_haplotype(genotype: str) -> int | None:
    is_phased_heterozygote = (
        "|" in genotype and Zygosity.of_genotype(genotype) is Zygosity.Heterozygous
    )
    if not is_phased_heterozygote:
        return None
    calls = alleles_of_genotype(genotype)
    return next(i for i, call in enumerate(calls) if call not in ("0", ".", ""))


def _allele_of_row(row: Row) -> Allele:
    return Allele(row["CHROM"], int(row["POS"]), row["REF"], row["ALT"])


def is_flagged_by_clinvar(row: Row) -> bool:
    # Also matches Conflicting_classifications_of_pathogenicity.
    return "athogenic" in row["CLNSIG"] or "athogenic" in row["CLNSIGCONF"]


def _clinvar_genes(row: Row) -> list[str]:
    """Symbols from a field such as 'BUB1B:701|BUB1B-PAK6:106821730'."""
    if _is_missing(row["GENEINFO"]):
        return []
    return [entry.split(":")[0] for entry in row["GENEINFO"].split("|")]


def _original_site(field: str) -> str | None:
    """Drop the allele index from bcftools' 'CHROM|POS|REF|ALTs|index'."""
    if _is_missing(field):
        return None
    return "|".join(field.split("|")[:4])


def _allele_depths(field: str) -> tuple[int, ...] | None:
    depths = field.split(",")
    if any(_is_missing(depth) for depth in depths):
        return None
    return tuple(map(int, depths))


def _largest_number(field: str) -> float | None:
    """Take the maximum across values separated by '&' or ','."""
    numbers = [float(part) for part in re.split("[&,]", field) if not _is_missing(part)]
    return max(numbers, default=None)


def _text(field: str) -> str:
    return "" if _is_missing(field) else field


def _is_missing(field: str) -> bool:
    return field in (".", "")
