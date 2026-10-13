import csv
from collections.abc import Sequence
from pathlib import Path

from mvapipe.ranking import Candidate, Evidence

CANDIDATE_COLUMNS = [
    "rank",
    "gene",
    "model",
    "score",
    "epcr",
    "genotype",
    "gene_prior",
    "phenotype",
    "disease",
    "phase",
]
VARIANT_COLUMNS = [
    "variant",
    "consequence",
    "zygosity",
    "clinvar",
    "clinvar_stars",
    "gnomad_af",
    "alphamissense",
    "impact",
    "filter",
    "gq",
    "dp",
]


def write_ranked_table(candidates: Sequence[Candidate], table_path: Path) -> None:
    header = CANDIDATE_COLUMNS + [
        f"v{number}_{column}" for number in (1, 2) for column in VARIANT_COLUMNS
    ]
    with table_path.open("w", encoding="utf-8", newline="") as table:
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for rank, candidate in enumerate(candidates, start=1):
            hypothesis = candidate.hypothesis
            phenotype = candidate.phenotype
            row: list[object] = [
                rank,
                hypothesis.gene,
                hypothesis.model.value,
                f"{candidate.score:.4f}",
                f"{candidate.epcr:.6g}",
                f"{hypothesis.genotype_evidence:.4f}",
                f"{hypothesis.gene_prior:.2f}",
                "" if phenotype is None else f"{phenotype.similarity:.3f}",
                ""
                if phenotype is None or phenotype.disease is None
                else phenotype.disease,
                hypothesis.phase.value,
            ]
            for evidence in sorted(hypothesis.variants, key=_position):
                row += _variant_cells(evidence)
            row += [""] * (len(header) - len(row))
            writer.writerow(row)


def _position(evidence: Evidence) -> int:
    return evidence.variant.allele.position


def _variant_cells(evidence: Evidence) -> list[object]:
    variant = evidence.variant
    return [
        variant.allele.to_string(),
        variant.consequence,
        evidence.zygosity.value,
        variant.clinvar_significance,
        evidence.review_stars,
        f"{variant.allele_frequency:.3g}",
        _formatted(variant.alphamissense, ".3f"),
        f"{evidence.impact:.2f}",
        variant.filter_status,
        _formatted(variant.genotype_quality, "g"),
        _formatted(variant.depth, "g"),
    ]


def _formatted(number: float | None, specification: str) -> str:
    return "" if number is None else format(number, specification)
