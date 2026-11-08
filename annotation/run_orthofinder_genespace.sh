#!/bin/bash

#SBATCH --job-name=orthofinder_gs
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=orthofinder_gs_%j.log

# Full all-vs-all OrthoFinder run on all 6 proteomes in orthofinder_input/,
# for GENESPACE. Kept separate from run_orthofinder.sh's --core/--assign run
# because GENESPACE needs the all-vs-all BLAST files for every genome pair.
# Conda env: orthofinder
#   built by annotation/setup_orthofinder.sh:
#   conda create -n orthofinder -c bioconda -c conda-forge orthofinder=3.1.5
# Usage: sbatch run_orthofinder_genespace.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

CONDA_BASE=/private/home/esteiger/miniconda3
ENV_NAME=orthofinder

source "${CONDA_BASE}/etc/profile.d/conda.sh"

# the aster package's activate.d script fails under set -u without this
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

conda activate "${ENV_NAME}"

INDIR=${BASE_DIR}/genome/annotation/orthofinder_input
OUTDIR=${BASE_DIR}/genome/annotation/orthofinder_results_genespace

# call the env's python explicitly; /bin precedes the env on PATH here
"${CONDA_PREFIX}/bin/python3" "${CONDA_PREFIX}/bin/orthofinder" \
    -f "${INDIR}" -o "${OUTDIR}" -t "${SLURM_CPUS_PER_TASK}" -a "${SLURM_CPUS_PER_TASK}"

echo "Done. Full all-vs-all results under: ${OUTDIR}/Results_*/"
echo "Next: point GENESPACE's genespace/orthofinder/ at this Results_*/ dir (see setup_genespace_run.sh)"
