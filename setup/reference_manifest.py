import argparse
import json
import re
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

from mvapipe.process import capture
from mvapipe.references import References
from mvapipe.vcf_steps import VEP_CACHE_VERSION

REFERENCE_DIRECTORIES = ["fasta", "clinvar", "alphamissense", "hpo", "gencc"]


def write_manifest(references: References) -> None:
    files = {
        path.relative_to(references.directory).as_posix(): {
            "bytes": path.stat().st_size,
        }
        for directory in REFERENCE_DIRECTORIES
        for path in sorted((references.directory / directory).iterdir())
    }
    clinvar_header = capture(["bcftools", "view", "-h", str(references.clinvar)])
    manifest = {
        "created": datetime.now(UTC).isoformat(timespec="seconds"),
        "files": files,
        "clinvar_file_date": _only_match(r"##fileDate=(\S+)", clinvar_header),
        "hpo_data_version": _only_match(
            r"data-version: (\S+)", references.hpo_obo.read_text(encoding="utf8")
        ),
        "vep_cache_version": VEP_CACHE_VERSION,
        "vep_plugins_revision": capture(
            ["git", "-C", str(references.vep_plugins), "rev-parse", "HEAD"]
        ).strip(),
        "tools": {
            "vep": _only_match(r"ensembl-vep\s*: (\S+)", capture(["vep", "--help"])),
            "bcftools": capture(["bcftools", "--version"]).splitlines()[0],
            "samtools": capture(["samtools", "--version"]).splitlines()[0],
        },
    }
    references.manifest.write_text(json.dumps(manifest, indent=1), encoding="utf8")


def _only_match(pattern: str, text: str) -> str:
    matches: list[str] = re.findall(pattern, text)
    (match,) = matches
    return match


class _Arguments(argparse.Namespace):
    refs: Path


def main(arguments: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="reference-manifest")
    parser.add_argument("--refs", type=Path, required=True)
    parsed = parser.parse_args(arguments, namespace=_Arguments())
    write_manifest(References(parsed.refs))


if __name__ == "__main__":
    main()
