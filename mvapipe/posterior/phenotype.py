import csv
import math
import re
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy import sparse

from mvapipe.ontology import Ontology, PhenotypeModel
from mvapipe.posterior.density_ratio import DensityRatio, fit_density_ratio

# Equal-width positive-similarity bins; missing annotations and zero have separate bins.
SIMILARITY_SCHEMES = {
    bins: tuple(round(index / bins, 4) for index in range(1, bins))
    for bins in (5, 10, 20, 40)
}

# Map (disease, current HPO term) to references from phenotype.hpoa.
AnnotationSources = dict[tuple[str, str], frozenset[str]]


class SimilarityIndex:
    """Compute `PhenotypeModel` similarities for all annotated genes in one pass."""

    def __init__(self, model: PhenotypeModel) -> None:
        self.ontology = model.ontology
        self.disease_terms = model.disease_terms
        genes_of_term: dict[str, set[str]] = defaultdict(set)
        closures: list[frozenset[str]] = []
        # Each pair identifies a row in the disease-ancestor matrix.
        self.pairs: list[tuple[str, str]] = []
        for gene, terms_of_disease in model.disease_terms.items():
            for disease, terms in terms_of_disease.items():
                closure = frozenset().union(*map(self.ontology.ancestors, terms))
                for ancestor in closure:
                    genes_of_term[ancestor].add(gene)
                closures.append(closure)
                self.pairs.append((gene, disease))
        gene_count = len(model.disease_terms)
        indexed_terms = sorted(genes_of_term)
        self._term_index = {term: index for index, term in enumerate(indexed_terms)}
        self._information_content = np.array(
            [-math.log(len(genes_of_term[term]) / gene_count) for term in indexed_terms]
        )
        self.genes = sorted(model.disease_terms)
        self.gene_index = {gene: index for index, gene in enumerate(self.genes)}
        self.pair_gene = np.array([self.gene_index[gene] for gene, _ in self.pairs])
        rows = [index for index, closure in enumerate(closures) for _ in closure]
        columns = [self._term_index[term] for closure in closures for term in closure]
        self._closures = sparse.csc_matrix(
            (np.ones(len(rows), dtype=bool), (rows, columns)),
            shape=(len(closures), len(indexed_terms)),
        )

    def score_pairs(self, patient_terms: Sequence[str]) -> NDArray[np.float64]:
        """Score nonempty patient terms against diseases in `pairs` order."""
        total = np.zeros(self._closures.shape[0])
        for term in patient_terms:
            ancestors = [
                self._term_index[ancestor]
                for ancestor in self.ontology.ancestors(term)
                if ancestor in self._term_index
            ]
            if not ancestors:
                continue
            values = self._information_content[ancestors]
            information_content = values.max()
            if information_content <= 0:
                continue
            shared = self._closures[:, ancestors].toarray() * values
            total += np.minimum(shared.max(axis=1) / information_content, 1.0)
        return total / len(patient_terms)

    def score_closure(
        self, closure: frozenset[str], patient_terms: Sequence[str]
    ) -> float:
        """Score nonempty patient terms against a disease's terms and ancestors."""
        total = 0.0
        for term in patient_terms:
            ancestors = [
                ancestor
                for ancestor in self.ontology.ancestors(term)
                if ancestor in self._term_index
            ]
            values = [self._information_content[self._term_index[a]] for a in ancestors]
            if not values or max(values) <= 0:
                continue
            shared = [
                self._information_content[self._term_index[ancestor]]
                for ancestor in ancestors
                if ancestor in closure
            ]
            total += min(max(shared, default=0.0) / max(values), 1.0)
        return total / len(patient_terms)

    def compute_gene_maxima(
        self, pair_scores: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Take each gene's maximum disease score, floored at 0, in `genes` order."""
        scores = np.zeros(len(self.genes))
        np.maximum.at(scores, self.pair_gene, pair_scores)
        return scores

    def score_genes(self, patient_terms: Sequence[str]) -> NDArray[np.float64]:
        """Score nonempty patient terms against genes in `genes` order."""
        return self.compute_gene_maxima(self.score_pairs(patient_terms))


def normalise_reference(reference: str) -> str:
    """Strip surrounding whitespace and standardise PubMed IDs as 'PMID:<number>'."""
    reference = reference.strip()
    match = re.fullmatch(r"PMID:\s*(\d+)", reference, flags=re.IGNORECASE)
    return f"PMID:{match.group(1)}" if match else reference


def find_pubmed_ids(references: Iterable[str]) -> frozenset[str]:
    """Extract unique PubMed references in 'PMID:<number>' form."""
    normalised = map(normalise_reference, references)
    return frozenset(
        reference for reference in normalised if reference.startswith("PMID:")
    )


def read_annotation_sources(hpoa_path: Path, ontology: Ontology) -> AnnotationSources:
    """Group phenotype.hpoa references by disease and current HPO term.

    Skip rows with the NOT qualifier and terms that are unknown or obsolete.
    """
    sources: dict[tuple[str, str], set[str]] = defaultdict(set)
    with hpoa_path.open(encoding="utf8", newline="") as lines:
        rows = csv.DictReader(
            (line for line in lines if not line.startswith("#")), delimiter="\t"
        )
        for row in rows:
            if row["qualifier"] == "NOT":
                continue
            term = ontology.resolve(row["hpo_id"])
            if term is not None:
                sources[row["database_id"], term] |= {
                    normalise_reference(reference)
                    for reference in row["reference"].split(";")
                    if reference.strip()
                }
    return {key: frozenset(references) for key, references in sources.items()}


class LeaveSourceOut:
    """Score genes after excluding annotations sourced only from the case being scored."""

    def __init__(self, similarity: SimilarityIndex, sources: AnnotationSources) -> None:
        self._similarity = similarity
        self._sources = sources
        self._terms_of_pair = [
            frozenset(similarity.disease_terms[gene][disease])
            for gene, disease in similarity.pairs
        ]
        self._pairs_of_publication: dict[str, list[int]] = defaultdict(list)
        for index, (_, disease) in enumerate(similarity.pairs):
            cited: set[str] = set()
            for term in self._terms_of_pair[index]:
                cited |= sources.get((disease, term), frozenset())
            for reference in cited:
                if reference.startswith("PMID:"):
                    self._pairs_of_publication[reference].append(index)

    def _find_kept_terms(
        self, index: int, publications: frozenset[str]
    ) -> frozenset[str]:
        disease = self._similarity.pairs[index][1]
        kept = set()
        for term in self._terms_of_pair[index]:
            references = self._sources.get((disease, term))
            if references is None or not references <= publications:
                kept.add(term)
        return frozenset(kept)

    def score_genes(
        self, patient_terms: Sequence[str], publications: frozenset[str]
    ) -> tuple[NDArray[np.float64], NDArray[np.bool_]]:
        """Return gene similarities and a mask of genes with retained annotations.

        Exclude terms whose source references are all in the case's `publications`.
        Return arrays in the similarity index's gene order.
        """
        similarity = self._similarity
        pair_scores = similarity.score_pairs(patient_terms)
        is_empty = np.zeros(len(pair_scores), dtype=bool)
        affected = {
            index
            for publication in publications
            for index in self._pairs_of_publication.get(publication, ())
        }
        for index in affected:
            kept = self._find_kept_terms(index, publications)
            if kept == self._terms_of_pair[index]:
                continue
            if not kept:
                pair_scores[index] = 0.0
                is_empty[index] = True
                continue
            closure = frozenset().union(*map(similarity.ontology.ancestors, kept))
            pair_scores[index] = similarity.score_closure(closure, patient_terms)
        is_annotated = np.zeros(len(similarity.genes), dtype=bool)
        np.logical_or.at(is_annotated, similarity.pair_gene, ~is_empty)
        return similarity.compute_gene_maxima(pair_scores), is_annotated


def build_similarity_features(
    similarity: NDArray[np.float64],
    is_annotated: NDArray[np.bool_],
    bounds: Sequence[float],
) -> NDArray[np.float64]:
    """Encode annotation absence, zero similarity and positive-similarity bins.

    Return one feature row per input, with columns in that order.
    """
    features = np.zeros((len(similarity), 2 + len(bounds) + 1))
    features[~is_annotated, 0] = 1
    features[is_annotated & (similarity <= 0), 1] = 1
    bins = np.searchsorted(bounds, similarity, side="left")
    is_positive = is_annotated & (similarity > 0)
    features[np.flatnonzero(is_positive), 2 + bins[is_positive]] = 1
    return features


@dataclass(frozen=True)
class PhenotypeTable:
    """Training features grouped by case, with its causal gene before sampled genes."""

    similarity: NDArray[np.float64]
    is_annotated: NDArray[np.bool_]
    is_causal: NDArray[np.bool_]
    fold: NDArray[np.int_]
    # Position within each case: 0 for the causal gene, then 1, 2, ... for sampled genes.
    order: NDArray[np.int_]


@dataclass(frozen=True)
class PhenotypeRatio:
    ratio: DensityRatio
    bounds: tuple[float, ...]

    def log_ratio(
        self, similarity: NDArray[np.float64], is_annotated: NDArray[np.bool_]
    ) -> NDArray[np.float64]:
        features = build_similarity_features(similarity, is_annotated, self.bounds)
        return self.ratio.log_ratio(features)


def fit_phenotype_ratio(
    table: PhenotypeTable, *, similarity_bins: int, regularisation: float
) -> PhenotypeRatio:
    """Fit a phenotype log density ratio between causal and sampled genes."""
    bounds = SIMILARITY_SCHEMES[similarity_bins]
    features = build_similarity_features(table.similarity, table.is_annotated, bounds)
    ratio = fit_density_ratio(
        positive=features[table.is_causal],
        negative=features[~table.is_causal],
        regularisation=regularisation,
    )
    return PhenotypeRatio(ratio, bounds)


@dataclass(frozen=True)
class GenePhenotype:
    """Phenotype log likelihood ratios for one patient, with an unannotated fallback."""

    log_ratio_of_gene: dict[str, float]
    # Fallback for genes without annotated diseases.
    unannotated_log_ratio: float

    def log_ratio(self, gene: str) -> float:
        return self.log_ratio_of_gene.get(gene, self.unannotated_log_ratio)


def score_patient_phenotype(
    ratio: PhenotypeRatio,
    similarity: SimilarityIndex,
    patient_terms: Sequence[str],
    *,
    leave_source_out: LeaveSourceOut | None = None,
    publications: frozenset[str] = frozenset(),
) -> GenePhenotype:
    """Compute phenotype log ratios in the similarity index's gene order.

    Supply `leave_source_out` and the case's PubMed IDs for literature evaluation.
    Ordinary patient ranking uses the full disease annotations.
    """
    if publications and leave_source_out is None:
        message = "publication exclusion requires annotation source references"
        raise ValueError(message)
    if leave_source_out is None:
        scores = similarity.score_genes(patient_terms)
        is_annotated = np.ones(len(scores), dtype=bool)
    else:
        scores, is_annotated = leave_source_out.score_genes(patient_terms, publications)
    log_ratios = ratio.log_ratio(scores, is_annotated)
    unannotated = ratio.log_ratio(np.zeros(1), np.zeros(1, dtype=bool))
    return GenePhenotype(
        log_ratio_of_gene=dict(zip(similarity.genes, log_ratios.tolist(), strict=True)),
        unannotated_log_ratio=float(unannotated[0]),
    )
