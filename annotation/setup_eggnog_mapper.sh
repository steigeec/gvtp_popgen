#!/bin/bash
# One-time setup. Builds conda env eggnog-mapper, plus its database:
#   conda create -n eggnog-mapper -c bioconda -c conda-forge eggnog-mapper=2.1.15
# Safe to re-run. Run on a login node (download only, no heavy compute).
#
# Usage: bash setup_eggnog_mapper.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

CONDA_BASE=/private/home/esteiger/miniconda3
ENV_NAME=eggnog-mapper
DATA_DIR=${BASE_DIR}/software/eggnog-mapper-data

source "${CONDA_BASE}/etc/profile.d/conda.sh"

if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo "Conda env '${ENV_NAME}' already exists — skipping creation."
else
    echo "Creating conda env '${ENV_NAME}'..."
    conda create -y -n "${ENV_NAME}" -c bioconda -c conda-forge eggnog-mapper
fi

conda activate "${ENV_NAME}"

# download_eggnog_data.py (2.1.15) points at eggnogdb.embl.de, which no
# longer resolves; the data is served from eggnog5.embl.de
DOWNLOAD_SCRIPT="$(which download_eggnog_data.py)"
if grep -q "eggnogdb.embl.de" "${DOWNLOAD_SCRIPT}"; then
    echo "Patching stale hostname (eggnogdb.embl.de -> eggnog5.embl.de) in ${DOWNLOAD_SCRIPT}"
    sed -i 's/eggnogdb\.embl\.de/eggnog5.embl.de/g' "${DOWNLOAD_SCRIPT}"
fi

mkdir -p "${DATA_DIR}"

# remove zero-byte files from failed downloads
find "${DATA_DIR}" -maxdepth 1 -type f -empty -delete

# core + taxonomic databases, tens of GB
if [ -f "${DATA_DIR}/eggnog.db" ] && [ -f "${DATA_DIR}/eggnog_proteins.dmnd" ]; then
    echo "eggNOG database already present at ${DATA_DIR} — skipping download."
else
    echo "Downloading eggNOG reference database to ${DATA_DIR}..."
    download_eggnog_data.py -y --data_dir "${DATA_DIR}"
fi

echo "Setup complete. Env: ${ENV_NAME}, data dir: ${DATA_DIR}"
echo "Next: sbatch run_eggnog_mapper.sh"
