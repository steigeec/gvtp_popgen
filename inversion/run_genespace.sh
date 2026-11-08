#!/bin/bash

#SBATCH --job-name=genespace
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --output=genespace_%j.log

# Run the GENESPACE pipeline (run_genespace.R): synteny blocks (MCScanX_h)
# across all 21 pairs (15 cross + 6 self) of the 6 staged species, using the
# precomputed OrthoFinder run at genespace/orthofinder.
# Outputs: genespace/{results,riparian,dotplots,pangenes}/
# Prerequisites: setup_genespace.sh (conda env, MCScanX, GENESPACE R package)
# and run_orthofinder_genespace.sh (all-vs-all OrthoFinder run).
# Conda env: genespace
#   built by annotation/setup_genespace.sh: conda env with r-base=4.3,
#   MCScanX and GENESPACE's Bioconductor/CRAN dependencies, then
#   devtools::install_github("jtlovell/GENESPACE") (GENESPACE 1.3.1)
#
# Usage: sbatch run_genespace.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

CONDA_BASE=/private/home/esteiger/miniconda3
ENV_NAME=genespace

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

cd ${BASE_DIR}/inversions/synteny
"${CONDA_PREFIX}/bin/Rscript" ${BASE_DIR}/gvtp_popgen/inversion/run_genespace.R
