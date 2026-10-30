import csv
import math
import re
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

PHENOTYPIC_ABNORMALITY = "HP:0000118"
PAST_MEDICAL_HISTORY = "HP:0032443"
FEATURE_BRANCHES = frozenset({PHENOTYPIC_ABNORMALITY, PAST_MEDICAL_HISTORY})

TERM_ID_PATTERN = re.compile(r"HP:\d{7}")


class Ontology:
    """HPO terms, parent relationships and alternate IDs loaded from an OBO file."""

    def __init__(self, obo_path: Path) -> None:
        self.term_names: dict[str, str] = {}
        self.parents: dict[str, set[str]] = defaultdict(set)
        self._primary_ids: dict[str, str] = {}
        self._obsolete_terms: set[str] = set()
        self._ancestor_cache: dict[str, frozenset[str]] = {}
        self._parse(obo_path)

    def _parse(self, obo_path: Path) -> None:
        term: str | None = None
        is_term_stanza = False
        with obo_path.open(encoding="utf8") as obo:
            for raw_line in obo:
                line = raw_line.rstrip("\n")
                if line.startswith("["):
                    term = None
                    is_term_stanza = line == "[Term]"
                elif not is_term_stanza:
                    continue
                elif line.startswith("id: "):
                    term = line.removeprefix("id: ").strip()
                elif term is None:
                    continue
                elif line.startswith("name: "):
                    self.term_names[term] = line.removeprefix("name: ").strip()
                elif line.startswith("is_a: "):
                    parent = line.removeprefix("is_a: ").split("!")[0].strip()
                    self.parents[term].add(parent)
                elif line.startswith("alt_id: "):
                    self._primary_ids[line.removeprefix("alt_id: ").strip()] = term
                elif line.startswith("is_obsolete: true"):
                    self._obsolete_terms.add(term)

    def resolve(self, term: str) -> str | None:
        """Resolve a primary or alternate HPO ID to its current primary ID.

        Return None when the term is unknown or obsolete.
        """
        primary = self._primary_ids.get(term, term)
        if primary in self.term_names and primary not in self._obsolete_terms:
            return primary
        return None

    def ancestors(self, term: str) -> frozenset[str]:
        """Return the term and all its ancestors."""
        cached = self._ancestor_cache.get(term)
        if cached is not None:
            return cached
        found = {term}
        pending = [term]
        while pending:
            for parent in self.parents.get(pending.pop(), ()):
                if parent not in found:
                    found.add(parent)
                    pending.append(parent)
        self._ancestor_cache[term] = frozenset(found)
        return self._ancestor_cache[term]

    def is_feature(self, term: str) -> bool:
        """Check whether a term describes a patient feature.

        Accept descendants of Phenotypic abnormality or Past medical history, excluding
        the two branch roots.
        """
        return term not in FEATURE_BRANCHES and not FEATURE_BRANCHES.isdisjoint(
            self.ancestors(term)
        )


@dataclass(frozen=True)
class PhenotypeMatch:
    similarity: float
    disease: str | None


NO_MATCH = PhenotypeMatch(similarity=0.0, disease=None)


class PhenotypeModel:
    """Score genes by mean normalised Resnik similarity to their diseases."""

    def __init__(self, *, obo_path: Path, genes_to_phenotype_path: Path) -> None:
        self.ontology = Ontology(obo_path)
        # Map gene symbols to disease IDs, then to current HPO terms.
        self.disease_terms = self._load_annotations(genes_to_phenotype_path)
        self._information_content = self._compute_information_content()
        self._closure_cache: dict[tuple[str, str], frozenset[str]] = {}

    def _load_annotations(
        self, genes_to_phenotype_path: Path
    ) -> dict[str, dict[str, set[str]]]:
        disease_terms: defaultdict[str, defaultdict[str, set[str]]] = defaultdict(
            lambda: defaultdict(set)
        )
        with genes_to_phenotype_path.open(encoding="utf8", newline="") as annotations:
            for row in csv.DictReader(annotations, delimiter="\t"):
                term = self.ontology.resolve(row["hpo_id"])
                if term is not None:
                    disease_terms[row["gene_symbol"]][row["disease_id"]].add(term)
        return {gene: dict(diseases) for gene, diseases in disease_terms.items()}

    def _compute_information_content(self) -> dict[str, float]:
        """Compute -log of the gene frequency of each term, including descendants."""
        genes_of_term: dict[str, set[str]] = defaultdict(set)
        for gene, terms_of_disease in self.disease_terms.items():
            for terms in terms_of_disease.values():
                for term in terms:
                    for ancestor in self.ontology.ancestors(term):
                        genes_of_term[ancestor].add(gene)
        gene_count = len(self.disease_terms)
        return {
            term: -math.log(len(genes) / gene_count)
            for term, genes in genes_of_term.items()
        }

    def _compute_closure(self, gene: str, disease: str) -> frozenset[str]:
        cached = self._closure_cache.get((gene, disease))
        if cached is not None:
            return cached
        closure: set[str] = set()
        for term in self.disease_terms[gene][disease]:
            closure |= self.ontology.ancestors(term)
        self._closure_cache[gene, disease] = frozenset(closure)
        return self._closure_cache[gene, disease]

    def _find_maximum_information_content(self, terms: frozenset[str]) -> float:
        return max(
            (self._information_content.get(term, 0.0) for term in terms), default=0.0
        )

    def score_gene(self, gene: str, patient_terms: Sequence[str]) -> PhenotypeMatch:
        """Find the gene's disease with the highest similarity to the patient terms.

        Return NO_MATCH if the gene has no annotated disease, the patient terms are empty
        or every similarity is zero.
        """
        if gene not in self.disease_terms or not patient_terms:
            return NO_MATCH
        patient = [
            (ancestors, self._find_maximum_information_content(ancestors))
            for ancestors in map(self.ontology.ancestors, patient_terms)
        ]
        best = NO_MATCH
        for disease in self.disease_terms[gene]:
            closure = self._compute_closure(gene, disease)
            total = sum(
                min(
                    self._find_maximum_information_content(ancestors & closure)
                    / information,
                    1.0,
                )
                for ancestors, information in patient
                if information > 0
            )
            similarity = total / len(patient)
            if similarity > best.similarity:
                best = PhenotypeMatch(similarity=similarity, disease=disease)
        return best


def read_term_file(term_path: Path, ontology: Ontology) -> list[str]:
    """Read unique, current HPO patient features in file order.

    Expect one HPO ID per line; ignore trailing text, '#' comments and blank lines.
    Raise ValueError for an invalid feature or a file with no terms.
    """
    terms: list[str] = []
    with term_path.open(encoding="utf8") as lines:
        for line_number, raw_line in enumerate(lines, start=1):
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            match = TERM_ID_PATTERN.match(line)
            term = ontology.resolve(match.group()) if match else None
            if term is None or not ontology.is_feature(term):
                message = (
                    f"{term_path}, line {line_number}: not a current HPO patient feature"
                )
                raise ValueError(message)
            if term not in terms:
                terms.append(term)
    if not terms:
        message = f"{term_path} has no terms"
        raise ValueError(message)
    return terms
