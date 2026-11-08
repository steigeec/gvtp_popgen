#!/bin/bash

#SBATCH --job-name=orthofinder
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=orthofinder_%j.log

# Add Helianthus anomalus to the existing 5-species orthogroups (Results_Jul15:
# Deinandra, Arabidopsis, H. annuus, Lactuca sativa, Solanum lycopersicum)
# with OrthoFinder 3's --core/--assign, without redoing the all-vs-all search.
# Output goes to a new Results_<date>/ next to Results_Jul15 (-o is not
# accepted with --core/--assign).
# orthofinder_input_assign/ must contain only the new species' proteome
# (prep_ncbi_datasets_proteome.sh GCA_047456115.1 Helianthus_anomalus).
# Requires: setup_orthofinder.sh (with numpy 2 patch).
# Conda env: orthofinder
#   built by annotation/setup_orthofinder.sh:
#   conda create -n orthofinder -c bioconda -c conda-forge orthofinder=3.1.5
# famsa crashes with "Illegal instruction" on the login node; run via sbatch.
# Usage: sbatch run_orthofinder.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

CONDA_BASE=/private/home/esteiger/miniconda3
ENV_NAME=orthofinder

source "${CONDA_BASE}/etc/profile.d/conda.sh"

# the aster package's activate.d script fails under set -u without this
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

conda activate "${ENV_NAME}"

CORE_RESULTS=${BASE_DIR}/genome/annotation/orthofinder_results/Results_Jul15
ASSIGN_DIR=${BASE_DIR}/genome/annotation/orthofinder_input_assign

# call the env's python explicitly; /bin precedes the env on PATH here
"${CONDA_PREFIX}/bin/python3" "${CONDA_PREFIX}/bin/orthofinder" \
    --core "${CORE_RESULTS}" --assign "${ASSIGN_DIR}" \
    -t "${SLURM_CPUS_PER_TASK}" -a "${SLURM_CPUS_PER_TASK}"

echo "Done. New Results_<date>/ created alongside ${CORE_RESULTS} under orthofinder_results/."
echo "Next: python3 build_tair_ortholog_go.py (once pointed at the new results dir)"
