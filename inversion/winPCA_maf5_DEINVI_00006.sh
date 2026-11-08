#!/bin/bash

#SBATCH --job-name=winPCA_maf5_00006
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --mem=64gb
#SBATCH --time=72:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/inversions/all_localities/winPCA_maf5_00006_%j.log

# 5%-MAF, GL-based windowed PCA of DEINVI_00006 (all localities), a sign-flip-
# corrected copy (<prefix>_flipped), and a PC1 plot coloured by fog.
# Usage: sbatch winPCA_maf5_DEINVI_00006.sh
# Conda env: winPCA
#   conda create -n winPCA -c conda-forge python=3.10 numpy=1.23.5 pandas \
#     numba scikit-allel plotly matplotlib scipy
#   pip install pcangsd   # needed for genotype-likelihood (GL) input
#   WinPCA itself: git clone https://github.com/MoritzBlumer/winpca
#   (run as software/winpca/winpca)

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/inversions/all_localities

BEAGLE="${BASE_DIR}/variant_calling/angsd/postParalog_4x.beagle.gz"
WIN="${BASE_DIR}/software/winpca/winpca"
CONDA_SH=/private/home/esteiger/miniconda3/etc/profile.d/conda.sh
META="${BASE_DIR}/inversions/metadata.tsv"
PLOT_SCRIPT="${BASE_DIR}/gvtp_popgen/inversion/small_chromplot.py"

CHROM='DEINVI_00006'
LENGTH=70516456
WINDOW=1000000
HALF=$((WINDOW/2))
START=${HALF}
END=$((LENGTH - HALF))

source ${CONDA_SH}
conda activate winPCA

PREFIX=${CHROM}_GL_maf5

# 5%-MAF-filtered, GL-based windowed PCA
if [ ! -f ${PREFIX}.stat.tsv.gz ]; then
    ${WIN} pca -v GL -t 40 -w ${WINDOW} --increment 10000 -m 0.05 \
        ${BEAGLE} ${CHROM}:1-${LENGTH} ${PREFIX}
fi

# Fix auto-polarization sign flips in PC1 on a copy (`winpca flip` overwrites
# its input). Each region is bounded by two clean flips (sample loadings r ~ -1
# with the neighbouring window). Guide-sample polarization is not used: the guides
# (inversion homozygotes) sit near 0 outside the inversion and it added flips
# elsewhere (84 vs 16 across the scaffold). The weak "flips" at 6.07/6.14 Mb
# (r = -0.57/-0.72) are PC1 changing identity, not sign, and are left alone.
FLIP_REGIONS='4070000-4080000,8500000,45830000-45910000,49210000-49830000,54330000-54370000,61260000-61290000'
FLIPPED=${PREFIX}_flipped
if [ ! -f ${FLIPPED}.stat.tsv.gz ]; then
    for f in ${PREFIX}.*.tsv.gz; do cp ${f} ${FLIPPED}${f#${PREFIX}}; done
    ${WIN} flip -c 1 -w ${FLIP_REGIONS} ${FLIPPED}
fi

# small-scale PC1 plot with a separate legend
python3 ${PLOT_SCRIPT} \
    ${PREFIX} ${CHROM}:${START}-${END} \
    -m ${META} -g fog --numeric \
    -o ${PREFIX}_fog.pc_1.png \
    --legend_out ${PREFIX}_fog.pc_1.legend.png
