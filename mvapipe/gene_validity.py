import csv
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum, IntEnum
from pathlib import Path


class Validity(IntEnum):
    NO_EVIDENCE = 0
    LIMITED = 1
    SUPPORTIVE = 2
    MODERATE = 3
    STRONG = 4
    DEFINITIVE = 5


class Inheritance(Enum):
    AUTOSOMAL_RECESSIVE = "AR"
    AUTOSOMAL_DOMINANT = "AD"
    X_LINKED = "XL"


CLASSIFICATION_VALIDITY = {
    "definitive": Validity.DEFINITIVE,
    "strong": Validity.STRONG,
    "moderate": Validity.MODERATE,
    "supportive": Validity.SUPPORTIVE,
    "limited": Validity.LIMITED,
    "disputed evidence": Validity.NO_EVIDENCE,
    "refuted evidence": Validity.NO_EVIDENCE,
    "no known disease relationship": Validity.NO_EVIDENCE,
    "animal model only": Validity.NO_EVIDENCE,
}


@dataclass(frozen=True)
class GeneValidity:
    """A gene's strongest GenCC classification overall and by inheritance mode."""

    best: Validity = Validity.NO_EVIDENCE
    by_inheritance: Mapping[Inheritance, Validity] = field(default_factory=dict)

    def for_inheritance(self, inheritance: Inheritance) -> Validity:
        return self.by_inheritance.get(inheritance, Validity.NO_EVIDENCE)


UNKNOWN_GENE = GeneValidity()


def parse_inheritance(title: str) -> Inheritance | None:
    """Parse a GenCC inheritance title, returning None for an unsupported mode."""
    title = title.lower()
    is_x_linked = "x-linked" in title
    if "recessive" in title and not is_x_linked:
        return Inheritance.AUTOSOMAL_RECESSIVE
    if "semidominant" in title or ("dominant" in title and not is_x_linked):
        return Inheritance.AUTOSOMAL_DOMINANT
    if is_x_linked:
        return Inheritance.X_LINKED
    return None


def load_gene_validity(gencc_path: Path) -> dict[str, GeneValidity]:
    """Load each gene's strongest classifications from a GenCC TSV export.

    Keep the strongest classification overall and for each supported inheritance mode.
    """
    best: dict[str, Validity] = {}
    by_inheritance: dict[str, dict[Inheritance, Validity]] = defaultdict(dict)
    with gencc_path.open(encoding="utf8", newline="") as submissions:
        for row in csv.DictReader(submissions, delimiter="\t"):
            symbol = row["gene_symbol"].strip()
            if not symbol:
                continue
            validity = CLASSIFICATION_VALIDITY[
                row["classification_title"].strip().lower()
            ]
            best[symbol] = max(best.get(symbol, Validity.NO_EVIDENCE), validity)
            inheritance = parse_inheritance(row["moi_title"])
            if inheritance is not None:
                classified = by_inheritance[symbol]
                previous = classified.get(inheritance, Validity.NO_EVIDENCE)
                classified[inheritance] = max(previous, validity)
    return {
        symbol: GeneValidity(best=validity, by_inheritance=by_inheritance.get(symbol, {}))
        for symbol, validity in best.items()
    }
