import csv
import gzip
import itertools
import re
from collections import defaultdict
from pathlib import Path

from mvapipe.frequency import FrequencyMode, frequency_of_row
from mvapipe.variant import Allele, Variant, Zygosity, split_genotype

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

# VEP impact levels, from most to least severe.
IMPACTS = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

Row = dict[str, str]


def load_variants(
    table_path: Path, *, frequency_mode: FrequencyMode = FrequencyMode.ALL
) -> list[Variant]:
    """Load variants with HIGH, MODERATE or LOW impact, or a ClinVar flag.

    Expect rows for each position to be adjacent. Keep one annotation per allele and gene,
    and use the highest allele frequency across all rows for that allele.
    """
    variants: list[Variant] = []
    seen: set[tuple[Allele, str]] = set()
    with gzip.open(table_path, "rt", encoding="utf8", newline="") as table:
        rows = csv.DictReader(table, delimiter="\t")
        if rows.fieldnames is None:
            message = f"{table_path} has no header"
            raise ValueError(message)
        for _, rows_of_position in itertools.groupby(
            rows, key=lambda row: (row["CHROM"], row["POS"])
        ):
            by_allele: dict[Allele, list[Row]] = defaultdict(list)
            for row in rows_of_position:
                by_allele[_parse_allele(row)].append(row)
            for allele_rows in by_allele.values():
                for variant in _build_variants(allele_rows, frequency_mode):
                    key = (variant.allele, variant.gene)
                    if key in seen:
                        message = (
                            "variant table repeats an allele outside its position group"
                        )
                        raise ValueError(message)
                    seen.add(key)
                    variants.append(variant)
    return variants


def _build_variants(rows: list[Row], frequency_mode: FrequencyMode) -> list[Variant]:
    call_columns = ("GT", "FILTER", "GQ", "DP", "AD", "PGT", "PID", "PS")
    calls = {tuple(row.get(column, ".") for column in call_columns) for row in rows}
    if len(calls) != 1:
        message = f"conflicting sample calls for {_parse_allele(rows[0]).to_string()}"
        raise ValueError(message)
    identifiers = dict.fromkeys(
        identifier
        for row in rows
        for identifier in row["ID"].split(";")
        if not _is_missing(identifier)
    )
    identifier = ";".join(identifiers) or "."
    is_flagged = is_flagged_by_clinvar(rows[0])
    clinvar_genes = _parse_clinvar_genes(rows[0]) if is_flagged else []
    # Restrict ClinVar evidence to its named genes if VEP annotates any of them.
    has_clinvar_gene = any(row["SYMBOL"] in clinvar_genes for row in rows)
    row_frequencies = [
        frequency
        for row in rows
        if (frequency := frequency_of_row(row, frequency_mode)) is not None
    ]
    allele_frequency = max(row_frequencies, default=None)

    by_gene: dict[str, Variant] = {}
    # Sort by severity so VEP's row order does not affect which annotation is kept.
    for row in sorted(rows, key=_severity_key):
        symbol = row["SYMBOL"]
        is_other_gene = (
            has_clinvar_gene and not _is_missing(symbol) and symbol not in clinvar_genes
        )
        is_flagged_here = is_flagged and not is_other_gene

        if not _is_missing(symbol):
            gene, is_gene_from_clinvar = symbol, False
        elif is_flagged_here and clinvar_genes:
            gene, is_gene_from_clinvar = clinvar_genes[0], True
        else:
            continue
        if row["IMPACT"] == "MODIFIER" and not is_flagged_here:
            continue

        # Prefer VEP's gene assignment over a gene inferred from ClinVar.
        existing = by_gene.get(gene)
        if existing is not None and (
            is_gene_from_clinvar or not existing.is_gene_from_clinvar
        ):
            continue
        by_gene[gene] = _parse_variant(
            row,
            gene=gene,
            is_gene_from_clinvar=is_gene_from_clinvar,
            keep_clinvar=not is_other_gene,
            allele_frequency=allele_frequency,
            identifier=identifier,
        )
    return list(by_gene.values())


def _severity_key(row: Row) -> tuple[int, str, str, str, str]:
    return (
        IMPACTS.index(row["IMPACT"]),
        row["SYMBOL"],
        row["Consequence"],
        row["EXON"],
        row["am_pathogenicity"],
    )


def _parse_variant(
    row: Row,
    *,
    gene: str,
    is_gene_from_clinvar: bool,
    keep_clinvar: bool,
    allele_frequency: float | None,
    identifier: str,
) -> Variant:
    phase_set, phase_haplotype = _parse_phase(row)
    return Variant(
        allele=_parse_allele(row),
        gene=gene,
        consequence=row["Consequence"],
        vep_impact=row["IMPACT"],
        exon=_blank_if_missing(row["EXON"]),
        genotype=row["GT"],
        filter_status=row["FILTER"],
        genotype_quality=_parse_largest_number(row.get("GQ", ".")),
        depth=_parse_largest_number(row.get("DP", ".")),
        allele_depths=_parse_allele_depths(row.get("AD", ".")),
        phase_set=phase_set,
        phase_haplotype=phase_haplotype,
        original_site=_parse_original_site(row["ORIG"]),
        allele_frequency=allele_frequency,
        alphamissense=_parse_largest_number(row["am_pathogenicity"]),
        clinvar_significance=_blank_if_missing(row["CLNSIG"]) if keep_clinvar else "",
        clinvar_review_status=_blank_if_missing(row["CLNREVSTAT"])
        if keep_clinvar
        else "",
        clinvar_conflict=_blank_if_missing(row["CLNSIGCONF"]) if keep_clinvar else "",
        identifier=identifier,
        is_gene_from_clinvar=is_gene_from_clinvar,
    )


def _parse_phase(row: Row) -> tuple[str | None, int | None]:
    """Return the phase set and the alternate allele's haplotype."""
    genotype = row["GT"]
    haplotype = _find_alternate_haplotype(genotype)
    physical_genotype, physical_set = row.get("PGT", "."), row.get("PID", ".")
    if not _is_missing(physical_set) and not _is_missing(physical_genotype):
        haplotypes = physical_genotype.split("|")
        original = row["ORIG"].split("|")
        alternate_index = original[4] if len(original) == 5 else "1"
        haplotype = (
            haplotypes.index(alternate_index)
            if len(haplotypes) == 2 and haplotypes.count(alternate_index) == 1
            else None
        )
        return f"PID:{physical_set}", haplotype
    phase_set = row.get("PS", ".")
    if not _is_missing(phase_set) and "|" in genotype:
        return f"PS:{row['CHROM']}:{phase_set}", haplotype
    return None, haplotype


def _find_alternate_haplotype(genotype: str) -> int | None:
    is_phased_heterozygote = (
        "|" in genotype and Zygosity.from_genotype(genotype) is Zygosity.HETEROZYGOUS
    )
    if not is_phased_heterozygote:
        return None
    calls = split_genotype(genotype)
    return next(i for i, call in enumerate(calls) if call not in ("0", ".", ""))


def _parse_allele(row: Row) -> Allele:
    return Allele(row["CHROM"], int(row["POS"]), row["REF"], row["ALT"])


def is_flagged_by_clinvar(row: Row) -> bool:
    """Check CLNSIG and CLNSIGCONF for pathogenicity text."""
    # Also matches Conflicting_classifications_of_pathogenicity.
    return "athogenic" in row["CLNSIG"] or "athogenic" in row["CLNSIGCONF"]


def _parse_clinvar_genes(row: Row) -> list[str]:
    """Extract gene symbols from ClinVar's 'symbol:id|symbol:id' field."""
    if _is_missing(row["GENEINFO"]):
        return []
    return [entry.split(":")[0] for entry in row["GENEINFO"].split("|")]


def _parse_original_site(field: str) -> str | None:
    """Drop the allele index from bcftools' 'CHROM|POS|REF|ALTs|index'."""
    if _is_missing(field):
        return None
    return "|".join(field.split("|")[:4])


def _parse_allele_depths(field: str) -> tuple[int, ...] | None:
    depths = field.split(",")
    if any(_is_missing(depth) for depth in depths):
        return None
    return tuple(map(int, depths))


def _parse_largest_number(field: str) -> float | None:
    """Return the largest '&'- or ','-separated value, or None if all are missing."""
    numbers = [float(part) for part in re.split("[&,]", field) if not _is_missing(part)]
    return max(numbers, default=None)


def _blank_if_missing(field: str) -> str:
    return "" if _is_missing(field) else field


def _is_missing(field: str) -> bool:
    return field in (".", "")
