from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class References:
    directory: Path

    @property
    def fasta(self) -> Path:
        return self.directory / "fasta" / "GRCh38.primary_assembly.fa.gz"

    @property
    def vep_cache(self) -> Path:
        return self.directory / "vep_cache"

    @property
    def vep_plugins(self) -> Path:
        return self.directory / "vep_plugins"

    @property
    def alphamissense(self) -> Path:
        return self.directory / "alphamissense" / "AlphaMissense_hg38.tsv.gz"

    @property
    def clinvar(self) -> Path:
        return self.directory / "clinvar" / "clinvar.vcf.gz"

    @property
    def gencc(self) -> Path:
        return self.directory / "gencc" / "gencc-submissions.tsv"

    @property
    def hpo_obo(self) -> Path:
        return self.directory / "hpo" / "hp.obo"

    @property
    def genes_to_phenotype(self) -> Path:
        return self.directory / "hpo" / "genes_to_phenotype.txt"

    @property
    def manifest(self) -> Path:
        return self.directory / "manifest.json"
