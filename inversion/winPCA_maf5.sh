#!/bin/bash

#SBATCH --job-name=winPCA_maf5
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --mem=64gb
#SBATCH --time=168:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/inversions/all_localities/winPCA_maf5_%x_%j.log

# 5%-MAF, GL-based windowed PCA (1 Mb windows, 10 kb step) for one scaffold,
# with the same settings as winPCA_maf5_DEINVI_00006.sh / _Strauss.sh.
# PCA only (DEINVI_0006_stacked.py reads the .tsv.gz outputs); sign flips are
# fixed afterwards with `winpca flip` on a copy.
# Conda env: winPCA
#   conda create -n winPCA -c conda-forge python=3.10 numpy=1.23.5 pandas \
#     numba scikit-allel plotly matplotlib scipy
#   pip install pcangsd   # needed for genotype-likelihood (GL) input
#   WinPCA itself: git clone https://github.com/MoritzBlumer/winpca
#   (run as software/winpca/winpca)
#
# Usage:
#   sbatch -J <name> winPCA_maf5.sh <scaffold> [<sample list> <tag>]
# e.g.
#   sbatch -J maf5_00017 winPCA_maf5.sh DEINVI_00017
#   sbatch -J maf5_00017_Strauss winPCA_maf5.sh DEINVI_00017 \
#       BASE_DIR/inversions/all_localities/strauss_samples.txt _Strauss
# Output prefix: <scaffold>_GL_maf5<tag>

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

CHROM=$1
SAMPLES=${2:-}
TAG=${3:-}

cd ${BASE_DIR}/inversions/all_localities

BEAGLE="${BASE_DIR}/variant_calling/angsd/postParalog_4x.beagle.gz"
WIN="${BASE_DIR}/software/winpca/winpca"
CONDA_SH=/private/home/esteiger/miniconda3/etc/profile.d/conda.sh
FAI="${BASE_DIR}/inversions/filtered_scaffolds.fai"
WINDOW=1000000

LENGTH=$(awk -v c="${CHROM}" '$1 == c {print $2}' ${FAI})
[[ -n "${LENGTH}" ]] || { echo "ERROR: ${CHROM} not in ${FAI}" >&2; exit 1; }

source ${CONDA_SH}
conda activate winPCA

PREFIX=${CHROM}_GL_maf5${TAG}
SUBSET_ARG=()
[[ -n "${SAMPLES}" ]] && SUBSET_ARG=(-s "${SAMPLES}")

echo "$(date): ${PREFIX}: ${CHROM}:1-${LENGTH}${SAMPLES:+, samples from ${SAMPLES}}"
if [ ! -f ${PREFIX}.stat.tsv.gz ]; then
    ${WIN} pca -v GL -t 40 -w ${WINDOW} --increment 10000 -m 0.05 "${SUBSET_ARG[@]}" \
        ${BEAGLE} ${CHROM}:1-${LENGTH} ${PREFIX}
else
    echo "${PREFIX}.stat.tsv.gz already exists -- skipping"
fi
echo "$(date): done"
