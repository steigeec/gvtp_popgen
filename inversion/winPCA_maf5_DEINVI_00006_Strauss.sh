#!/bin/bash

#SBATCH --job-name=winPCA_maf5_00006_Strauss
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --mem=64gb
#SBATCH --time=120:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/inversions/all_localities/winPCA_maf5_00006_Strauss_%j.log

# Same as winPCA_maf5_DEINVI_00006.sh (5%-MAF, GL-based, 1 Mb windows / 10 kb
# step), but Strauss samples only, so the stacked DEINVI_00006 figure can show
# every individual without Strauss (277 of 339 samples) swamping the other
# localities. strauss_samples.txt = beagle samples with Clusters == Strauss
# in metadata.tsv.

cd /private/groups/meyerlab/esteiger/inversions/all_localities

BEAGLE='/private/groups/meyerlab/esteiger/variant_calling/angsd/postParalog_4x.beagle.gz'
WIN='/private/groups/meyerlab/esteiger/software/winpca/winpca'
CONDA_SH=/private/home/esteiger/miniconda3/etc/profile.d/conda.sh
SAMPLES='/private/groups/meyerlab/esteiger/inversions/all_localities/strauss_samples.txt'

CHROM='DEINVI_00006'
LENGTH=70516456
WINDOW=1000000

source ${CONDA_SH}
conda activate winPCA

PREFIX=${CHROM}_GL_maf5_Strauss

if [ ! -f ${PREFIX}.stat.tsv.gz ]; then
    ${WIN} pca -v GL -t 40 -w ${WINDOW} --increment 10000 -m 0.05 -s ${SAMPLES} \
        ${BEAGLE} ${CHROM}:1-${LENGTH} ${PREFIX}
fi

# Fix auto-polarization sign flips in PC1 (on a copy -- `winpca flip`
# overwrites its input), as for winPCA_maf5_DEINVI_00006.sh. Regions were
# found from this run's output: every window whose sample loadings sit on the
# other side of a clean flip (r < -0.9 with the neighbouring window). None fall
# inside INV1 (43-49.3 Mb). Weak "flips" left alone (PC1 changing identity,
# not sign): 61.46 Mb (r = -0.88), 61.51 (-0.76), 63.16 (-0.06).
FLIP_REGIONS='21130000-37700000,37860000-37950000,61580000-62100000,62170000-63250000,64040000-64060000,64100000,67380000-67610000,67740000-70010000'
FLIPPED=${PREFIX}_flipped
if [ ! -f ${FLIPPED}.stat.tsv.gz ]; then
    for f in ${PREFIX}.*.tsv.gz; do cp ${f} ${FLIPPED}${f#${PREFIX}}; done
    ${WIN} flip -c 1 -w ${FLIP_REGIONS} ${FLIPPED}
fi
