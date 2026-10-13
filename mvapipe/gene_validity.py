import csv
from dataclasses import dataclass, field
from enum import Enum, IntEnum
from pathlib import Path


class Validity(IntEnum):
    No_evidence = 0
    Limited = 1
    Supportive = 2
    Moderate = 3
    Strong = 4
    Definitive = 5


class Inheritance(Enum):
    Autosomal_recessive = "AR"
    Autosomal_dominant = "AD"
    X_linked = "XL"


CLASSIFICATION_VALIDITY = {
    "definitive": Validity.Definitive,
    "strong": Validity.Strong,
    "moderate": Validity.Moderate,
    "supportive": Validity.Supportive,
    "limited": Validity.Limited,
    "disputed evidence": Validity.No_evidence,
    "refuted evidence": Validity.No_evidence,
    "no known disease relationship": Validity.No_evidence,
    "animal model only": Validity.No_evidence,
}


@dataclass
class Gene_validity:
    """The strongest GenCC gene-disease classification, overall and per inheritance."""

    best: Validity = Validity.No_evidence
    by_inheritance: dict[Inheritance, Validity] = field(default_factory=dict)

    def for_inheritance(self, inheritance: Inheritance) -> Validity:
        return self.by_inheritance.get(inheritance, Validity.No_evidence)


UNKNOWN_GENE = Gene_validity()


def inheritance_of_title(title: str) -> Inheritance | None:
    title = title.lower()
    is_x_linked = "x-linked" in title
    if "recessive" in title and not is_x_linked:
        return Inheritance.Autosomal_recessive
    if "semidominant" in title or ("dominant" in title and not is_x_linked):
        return Inheritance.Autosomal_dominant
    if is_x_linked:
        return Inheritance.X_linked
    return None


def load_gene_validity(gencc_path: Path) -> dict[str, Gene_validity]:
    csv.field_size_limit(1 << 30)
    genes: dict[str, Gene_validity] = {}
    with gencc_path.open(encoding="utf8", newline="") as submissions:
        for row in csv.DictReader(submissions, delimiter="\t"):
            symbol = row["gene_symbol"].strip()
            if not symbol:
                continue
            validity = CLASSIFICATION_VALIDITY[
                row["classification_title"].strip().lower()
            ]
            gene = genes.setdefault(symbol, Gene_validity())
            gene.best = max(gene.best, validity)
            inheritance = inheritance_of_title(row["moi_title"])
            if inheritance is not None:
                previous = gene.for_inheritance(inheritance)
                gene.by_inheritance[inheritance] = max(previous, validity)
    return genes
