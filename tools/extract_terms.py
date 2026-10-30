import argparse
import re
import zipfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree

from mvapipe.ontology import Ontology
from mvapipe.references import References

WORD_NAMESPACE = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
CODE_PATTERN = re.compile(r"\bHP[:_ ]\s?(\d{7})\b", re.IGNORECASE)
ABSENCE_PATTERN = re.compile(
    r"\b(no|not|absent|absence|without|negative|excluded|ruled out|denies|denied)\b",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class Extraction:
    terms: list[str]
    # One-based positions in the supplied line sequence.
    lines_stating_absence: list[int]
    lines_with_unusable_code: list[int]


def lines_of_document(document_path: Path) -> list[str]:
    """Read nonblank text lines or DOCX paragraphs, including paragraphs in tables."""
    if document_path.suffix.lower() != ".docx":
        lines = document_path.read_text(encoding="utf-8-sig").splitlines()
    else:
        with zipfile.ZipFile(document_path) as archive:
            root = ElementTree.fromstring(archive.read("word/document.xml"))  # noqa: S314
        lines = [
            "".join(run.text or "" for run in paragraph.iter(f"{WORD_NAMESPACE}t"))
            for paragraph in root.iter(f"{WORD_NAMESPACE}p")
        ]
    return [line.strip() for line in lines if line.strip()]


def extract_terms(lines: Sequence[str], ontology: Ontology) -> Extraction:
    """Extract unique current HPO features, excluding lines with negation words.

    Remove HPO codes and term names before checking for negation. Report one-based input
    positions for negated lines and lines containing unusable codes.
    """
    extraction = Extraction([], [], [])
    for number, line in enumerate(lines, start=1):
        codes = [f"HP:{digits}" for digits in CODE_PATTERN.findall(line)]
        terms = [
            term
            for code in codes
            if (term := ontology.resolve(code)) is not None and ontology.is_feature(term)
        ]
        if len(terms) < len(codes):
            extraction.lines_with_unusable_code.append(number)
        if not terms:
            continue
        # Remove term names so features such as "Absent speech" do not imply negation.
        remainder = CODE_PATTERN.sub(" ", line)
        for term in terms:
            name = re.escape(ontology.term_names[term])
            remainder = re.sub(name, " ", remainder, flags=re.IGNORECASE)
        if ABSENCE_PATTERN.search(remainder):
            extraction.lines_stating_absence.append(number)
        else:
            extraction.terms.extend(
                term for term in terms if term not in extraction.terms
            )
    return extraction


def write_term_file(terms: Sequence[str], ontology: Ontology, term_path: Path) -> None:
    with term_path.open("w", encoding="utf8", newline="\n") as lines:
        for term in terms:
            lines.write(f"{term}  # {ontology.term_names[term]}\n")


class _Arguments(argparse.Namespace):
    document: Path
    refs: Path
    out: Path


def main(arguments: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="extract-terms")
    parser.add_argument("document", type=Path)
    parser.add_argument("--refs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parsed = parser.parse_args(arguments, namespace=_Arguments())
    ontology = Ontology(References(parsed.refs).hpo_obo)
    extraction = extract_terms(lines_of_document(parsed.document), ontology)
    write_term_file(extraction.terms, ontology, parsed.out)
    print(f"{len(extraction.terms)} terms written to {parsed.out}")
    print(
        f"lines not used because they state absence: {extraction.lines_stating_absence}"
    )
    print(
        "lines with a code that is not a current patient feature: "
        f"{extraction.lines_with_unusable_code}"
    )
    return 0 if extraction.terms else 1


if __name__ == "__main__":
    raise SystemExit(main())
