#!/bin/bash

#SBATCH --job-name=winPCA
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=40
#SBATCH --mem=64gb
#SBATCH --time=48:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/inversions/all_localities/winPCA_GT_%j.log

# GT-based windowed PCA and per-sample heterozygosity (hetp) per scaffold
# (1 Mb windows, 10 kb step), run with GNU parallel.
# Usage: sbatch [--export=ALL,SUBSET=...,TAG=...,SCAFFOLDS=...] winPCA_GT.sh
# Conda env: winPCA
#   conda create -n winPCA -c conda-forge python=3.10 numpy=1.23.5 pandas \
#     numba scikit-allel plotly matplotlib scipy
#   pip install pcangsd   # needed for genotype-likelihood (GL) input
#   WinPCA itself: git clone https://github.com/MoritzBlumer/winpca
#   (run as software/winpca/winpca)

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/inversions/all_localities

FAI="${BASE_DIR}/inversions/filtered_scaffolds.fai"
# Not MAF-filtered, so per-sample heterozygosity (hetp) isn't biased by
# dropping rare (mostly heterozygous) variants
VCF="${BASE_DIR}/variant_calling/angsd/postNGSparalog_beforeMAF_auto.vcf"
# Samples to include: >=8x coverage only (hard-called GTs undercount hets at
# lower coverage). Override at submit time, e.g. for Strauss only:
#   sbatch --export=ALL,SUBSET=BASE_DIR/inversions/all_localities/strauss_samples_8x.txt,TAG=_Strauss winPCA_GT.sh
# (strauss_samples_8x.txt = bamlist_8x.txt intersected with strauss_samples.txt)
SUBSET=${SUBSET:-"${BASE_DIR}/mapping/merged/bamlist_8x.txt"}
# appended to output prefixes/logs so runs on different sample sets don't
# collide (an existing log makes the skip check below skip the scaffold)
TAG=${TAG:-'_all8x'}
WIN="${BASE_DIR}/software/winpca/winpca"
CONDA_SH=/private/home/esteiger/miniconda3/etc/profile.d/conda.sh
META="${BASE_DIR}/inversions/metadata.tsv"
COLORS='Strauss:009e73,CDFW:e69f00,Dangermond:56b4e9,VAFB_Arguello:d55e00,VAFB_Lion_Head:f0e442,Point_Conception:0072b2'
# Scaffolds to run (regex on the FAI's first column); override at submit time,
# e.g. sbatch --export=ALL,SCAFFOLDS='^DEINVI_00017$' winPCA_GT.sh
SCAFFOLDS=${SCAFFOLDS:-'^DEINVI_00006$'}

export FAI
export COLORS
export CONDA_SH
export WIN
export VCF
export SUBSET
export META
export TAG

export PATH=/usr/bin:/bin:$PATH

awk -v re="${SCAFFOLDS}" '$1 ~ re' ${FAI} | parallel -j 5  --colsep '\t' --env COLORS --env META --env PATH --env CONDA_SH --env WIN --env VCF --env SUBSET --env TAG "

    WINDOW=1000000
    INCREMENT=10000
    HALF=\$((WINDOW/2))
    END=\$(({2} - HALF))
    source \${CONDA_SH}
    conda activate winPCA

    # winpca streams (and splits) every line of the VCF, twice, so first pull
    # this scaffold out of the full ~1.8TB VCF in one fast pass (stops reading
    # once past the scaffold, assuming the VCF is sorted by scaffold)
    SCAF_VCF={1}_beforeMAF.vcf.gz
    if [ ! -f \${SCAF_VCF} ]; then
        echo \"[parallel] extracting {1} from \${VCF}\"
        LC_ALL=C awk 'substr(\$0,1,1)==\"#\" {print; next}
                      index(\$0, \"{1}\t\")==1 {f=1; print; next}
                      f {exit}' \${VCF} | gzip -1 > \${SCAF_VCF}.tmp && mv \${SCAF_VCF}.tmp \${SCAF_VCF}
    fi
    LOG={1}_GT\${TAG}.log

    # --- SBATCH log diagnostics (stdout) ---
    echo \"[parallel] START scaffold={1} length={2} END=\${END}\"
    echo \"  stat file ({1}_GT\${TAG}.stat.tsv.gz): \$([ -f {1}_GT\${TAG}.stat.tsv.gz ] && echo EXISTS || echo MISSING)\"
    echo \"  log  file ({1}_GT\${TAG}.log):         \$([ -f {1}_GT\${TAG}.log ]         && echo EXISTS || echo MISSING)\"
    echo \"  png  file ({1}_GT\${TAG}_Clusters.pc_1.png): \$([ -f {1}_GT\${TAG}_Clusters.pc_1.png ] && echo EXISTS || echo MISSING)\"

    # skip if a log file already exists
    if [ -f {1}_GT\${TAG}.log ]; then
        echo 'skipping because log file already exists' >> {1}_GT\${TAG}.log 2>&1
        echo \"[parallel] SKIPPING {1} (log file already exists)\"
        exit 0
    fi

    # written after the skip check so it doesn't trigger the skip
    echo \"DEBUG: Processing scaffold {1} with length {2} END=\${END}\" >> {1}_GT\${TAG}.log 2>&1

    # GT-based PCA, >=8x samples (SUBSET)
    if [ ! -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        \${WIN} pca -v GT -t 40 -w 1000000 -s \${SUBSET} --increment 10000 \${SCAF_VCF} {1}:1-{2} {1}_GT\${TAG} >> {1}\${TAG}.log 2>&1
    fi

    if [ ! -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        echo 'Error: {1}_GT\${TAG}.stat.tsv.gz for genotype-based work not found!' >> {1}\${TAG}.log 2>&1
        exit 1
    fi

"
