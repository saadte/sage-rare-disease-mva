#!/usr/bin/env bash
# Usage: MVA_ROOT=DIRECTORY setup_env.sh
set -euo pipefail

REPOSITORY="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REFS="${MVA_ROOT:?}/refs"
ENVIRONMENT="${MVA_ROOT}/env"
MARKERS="${MVA_ROOT}/.setup"
MICROMAMBA="${MVA_ROOT}/bin/micromamba"
# shellcheck source=../config/resources.env
source "${REPOSITORY}/config/resources.env"
mkdir -p "${MVA_ROOT}/bin" "${REFS}" "${MARKERS}"
cd "${REPOSITORY}"

in_environment() { "${MICROMAMBA}" run -p "${ENVIRONMENT}" "$@"; }

step() {
  local name="$1"
  if [[ -e "${MARKERS}/${name}" ]]; then
    echo "[skip] ${name}"
    return
  fi
  echo "[run ] ${name}"
  "${name}"
  touch "${MARKERS}/${name}"
}

download() {
  mkdir -p "$(dirname "$2")"
  in_environment aria2c --quiet --continue --split=8 --max-connection-per-server=8 \
    --dir="$(dirname "$2")" --out="$(basename "$2")" "$1"
}

install_micromamba() {
  curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest |
    tar -xj -C "${MVA_ROOT}" bin/micromamba
}

create_environment() {
  rm -rf "${ENVIRONMENT}"
  "${MICROMAMBA}" create --yes --quiet -p "${ENVIRONMENT}" -f environment.lock.yml
}

download_vep_cache() {
  in_environment s5cmd --no-sign-request --log error --numworkers 64 \
    sync "${VEP_CACHE_S3}" "${REFS}/vep_cache/"
}

download_vep_plugins() {
  rm -rf "${REFS}/vep_plugins"
  in_environment git init --quiet "${REFS}/vep_plugins"
  in_environment git -C "${REFS}/vep_plugins" remote add origin "${VEP_PLUGINS_REPOSITORY}"
  in_environment git -C "${REFS}/vep_plugins" fetch --quiet --depth 1 origin "${VEP_PLUGINS_REVISION}"
  in_environment git -C "${REFS}/vep_plugins" checkout --quiet --detach FETCH_HEAD
}

download_fasta() {
  local gzipped="${REFS}/fasta/download.fa.gz"
  local fasta="${REFS}/fasta/GRCh38.primary_assembly.fa.gz"
  local threads
  threads="$(nproc)"
  download "${FASTA_URL}" "${gzipped}"
  # Recompress as bgzip so samtools can index the FASTA for random access.
  gzip -dc "${gzipped}" | in_environment bgzip -@ "${threads}" -c > "${fasta}"
  in_environment samtools faidx "${fasta}"
  rm "${gzipped}"
}

download_clinvar() {
  download "${CLINVAR_URL}" "${REFS}/clinvar/clinvar.vcf.gz"
  in_environment tabix -p vcf "${REFS}/clinvar/clinvar.vcf.gz"
}

download_alphamissense() {
  download "${ALPHAMISSENSE_URL}" "${REFS}/alphamissense/AlphaMissense_hg38.tsv.gz"
  # Index chromosome column 1 and position column 2, skipping the single header line.
  in_environment tabix -s 1 -b 2 -e 2 -S 1 "${REFS}/alphamissense/AlphaMissense_hg38.tsv.gz"
}

download_hpo() {
  download "${HPO_URL}/hp.obo" "${REFS}/hpo/hp.obo"
  download "${HPO_URL}/genes_to_phenotype.txt" "${REFS}/hpo/genes_to_phenotype.txt"
}

install_gencc() {
  mkdir -p "${REFS}/gencc"
  cp "${REPOSITORY}/config/gencc-2026-10-13.tsv" "${REFS}/gencc/gencc-submissions.tsv"
}

install_gene_folds() {
  mkdir -p "${REFS}/phenopackets"
  cp "${REPOSITORY}/config/gene_folds.tsv" "${REFS}/phenopackets/gene_folds.tsv"
}

download_phenotype_annotations() {
  download "${HPO_URL}/phenotype.hpoa" "${REFS}/hpo/phenotype.hpoa"
}

download_phenopackets() {
  download "${PHENOPACKETS_URL}" "${REFS}/phenopackets/all_phenopackets.zip"
}

download_gene_constraint() {
  download "${GNOMAD_CONSTRAINT_URL}" "${REFS}/gnomad/gnomad.v4.1.constraint_metrics.tsv"
}

prepare_training() {
  in_environment python -m mvapipe prepare-training --refs "${REFS}"
}

step install_micromamba
step create_environment
step download_vep_cache
step download_vep_plugins
step download_fasta
step download_clinvar
step download_alphamissense
step download_hpo
step install_gencc
step download_phenotype_annotations
step download_phenopackets
step download_gene_constraint
step install_gene_folds
step prepare_training
