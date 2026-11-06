#!/bin/bash

#SBATCH --job-name=winPCA
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=40
#SBATCH --mem=64gb
#SBATCH --time=48:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/inversions/all_localities/winPCA_GT_%j.log

cd /private/groups/meyerlab/esteiger/inversions/all_localities

FAI='/private/groups/meyerlab/esteiger/inversions/filtered_scaffolds.fai'
# Not MAF-filtered, so per-sample heterozygosity (hetp) isn't biased by
# dropping rare (mostly heterozygous) variants
VCF='/private/groups/meyerlab/esteiger/variant_calling/angsd/postNGSparalog_beforeMAF_auto.vcf'
# Samples to include: >=8x coverage only (hard-called GTs undercount hets at
# lower coverage). Override at submit time, e.g. for Strauss only:
#   sbatch --export=ALL,SUBSET=/private/groups/meyerlab/esteiger/inversions/all_localities/strauss_samples_8x.txt,TAG=_Strauss winPCA_GT.sh
# (strauss_samples_8x.txt = bamlist_8x.txt intersected with strauss_samples.txt)
SUBSET=${SUBSET:-'/private/groups/meyerlab/esteiger/mapping/merged/bamlist_8x.txt'}
# appended to output prefixes/logs so runs on different sample sets don't
# collide (an unsuffixed DEINVI_00006_GT.log from a cancelled run would make
# the skip check below skip the scaffold)
TAG=${TAG:-'_all8x'}
WIN='/private/groups/meyerlab/esteiger/software/winpca/winpca'
CONDA_SH=/private/home/esteiger/miniconda3/etc/profile.d/conda.sh
META='/private/groups/meyerlab/esteiger/inversions/metadata.tsv'
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

# Before the parallel command, set up the full PATH
export PATH=/usr/bin:/bin:$PATH


awk -v re="${SCAFFOLDS}" '$1 ~ re' ${FAI} | parallel -j 5  --colsep '\t' --env COLORS --env META --env PATH --env CONDA_SH --env WIN --env VCF --env SUBSET --env TAG "

    WINDOW=1000000
    INCREMENT=10000
    HALF=\$((WINDOW/2))
    END=\$(({2} - HALF))
    # Source conda setup in each subshell
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
    # Create log file for a given chrom
    LOG={1}_GT\${TAG}.log

    # --- SBATCH log diagnostics (stdout) ---
    echo \"[parallel] START scaffold={1} length={2} END=\${END}\"
    echo \"  stat file ({1}_GT\${TAG}.stat.tsv.gz): \$([ -f {1}_GT\${TAG}.stat.tsv.gz ] && echo EXISTS || echo MISSING)\"
    echo \"  log  file ({1}_GT\${TAG}.log):         \$([ -f {1}_GT\${TAG}.log ]         && echo EXISTS || echo MISSING)\"
    echo \"  png  file ({1}_GT\${TAG}_Clusters.pc_1.png): \$([ -f {1}_GT\${TAG}_Clusters.pc_1.png ] && echo EXISTS || echo MISSING)\"

    # Skip if there's already an out file
    if [ -f {1}_GT\${TAG}.log ]; then
        echo 'skipping because log file already exists' >> {1}_GT\${TAG}.log 2>&1
        echo \"[parallel] SKIPPING {1} (log file already exists)\"
        exit 0
    fi

    # Per-scaffold debug header (written after skip check so it doesn't trigger the skip)
    echo \"DEBUG: Processing scaffold {1} with length {2} END=\${END}\" >> {1}_GT\${TAG}.log 2>&1

    # GT-based PCA
    
    # Restrict to 8x+ only... 
    if [ ! -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        \${WIN} pca -v GT -t 40 -w 1000000 -s \${SUBSET} --increment 10000 \${SCAF_VCF} {1}:1-{2} {1}_GT\${TAG} >> {1}\${TAG}.log 2>&1
    fi

    # Ensure the stat file exists before using it
    if [ ! -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        echo 'Error: {1}_GT\${TAG}.stat.tsv.gz for genotype-based work not found!' >> {1}\${TAG}.log 2>&1
        exit 1
    fi

    # chromplot section commented out: the stacked figure (DEINVI_0006_stacked.py)
    # reads the .tsv.gz outputs directly, and these retries (up to 5 h each)
    # kept the job running for many hours after winpca had finished.

    # GRADUALLY TRY LARGER AND LARGER INTERVALS TO ACHIEVE SUCCESSFUL PLOTTING... 

    # PC plots
    # if [ ! -f {1}_GT\${TAG}_Clusters.pc_1.png ] && [ -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        # echo 'attempting for no interval' >> {1}_GT\${TAG}.log 2>&1
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.pc_1.png {1}_GT\${TAG}_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g Clusters -c \${COLORS} -p PC2 >> {1}_GT\${TAG}.log 2>&1
        #mv {1}_GT\${TAG}.pc_2.png {1}_GT\${TAG}_Clusters.pc_2.png 
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -p het >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.het.png {1}_GT\${TAG}_Clusters.het.png 
    # fi

    # if [ ! -f {1}_GT\${TAG}_Clusters.pc_1.png ] && [ -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        # echo 'attempting for interval of 5' >> {1}_GT\${TAG}.log 2>&1
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -i 5 >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.pc_1.png {1}_GT\${TAG}_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g Clusters -c \${COLORS} -p PC2 -i 5>> {1}_GT\${TAG}.log 2>&1
        #mv {1}_GT\${TAG}.pc_2.png {1}_GT\${TAG}_Clusters.pc_2.png 
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -p het -i 5 >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.het.png {1}_GT\${TAG}_Clusters.het.png 
    # fi

    # if [ ! -f {1}_GT\${TAG}_Clusters.pc_1.png ] && [ -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        # echo 'attempting for interval of 10' >> {1}_GT\${TAG}.log 2>&1
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -i 10 >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.pc_1.png {1}_GT\${TAG}_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g Clusters -c \${COLORS} -p PC2 -i 10 >> {1}_GT\${TAG}.log 2>&1
        #mv {1}_GT\${TAG}.pc_2.png {1}_GT\${TAG}_Clusters.pc_2.png 
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -p het -i 10 >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.het.png {1}_GT\${TAG}_Clusters.het.png 
    # fi

    # if [ ! -f {1}_GT\${TAG}_Clusters.pc_1.png ] && [ -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        # echo 'attempting for interval of 15' >> {1}_GT\${TAG}.log 2>&1
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -i 15 >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.pc_1.png {1}_GT\${TAG}_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g Clusters -c \${COLORS} -p PC2 -i 15>> {1}_GT\${TAG}.log 2>&1
        #mv {1}_GT\${TAG}.pc_2.png {1}_GT\${TAG}_Clusters.pc_2.png 
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -p het -i 15 >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.het.png {1}_GT\${TAG}_Clusters.het.png 
    # fi

    # if [ ! -f {1}_GT\${TAG}_Clusters.pc_1.png ] && [ -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        # echo 'attempting for interval of 20' >> {1}_GT\${TAG}.log 2>&1
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -i 20 >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.pc_1.png {1}_GT\${TAG}_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g Clusters -c \${COLORS} -p PC2 -i 20 >> {1}_GT\${TAG}.log 2>&1
        #mv {1}_GT\${TAG}.pc_2.png {1}_GT\${TAG}_Clusters.pc_2.png 
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -p het -i 20 >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.het.png {1}_GT\${TAG}_Clusters.het.png 
    # fi

    # if [ ! -f {1}_GT\${TAG}_Clusters.pc_1.png ] && [ -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        # echo 'attempting for interval of 25' >> {1}_GT\${TAG}.log 2>&1
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -i 25 >> {1}\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.pc_1.png {1}_GT\${TAG}_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g Clusters -c \${COLORS} -p PC2 -i 25>> {1}\${TAG}.log 2>&1
        #mv {1}_GT\${TAG}.pc_2.png {1}_GT\${TAG}_Clusters.pc_2.png 
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -p het -i 25 >> {1}\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.het.png {1}_GT\${TAG}_Clusters.het.png 
    # fi

    # if [ ! -f {1}_GT\${TAG}_Clusters.pc_1.png ] && [ -f {1}_GT\${TAG}.stat.tsv.gz ]; then
        # echo 'attempting for interval of 30' >> {1}_GT\${TAG}.log 2>&1
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -i 30 >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.pc_1.png {1}_GT\${TAG}_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g Clusters -c \${COLORS} -p PC2 -i 30>> {1}_GT\${TAG}.log 2>&1
        #mv {1}_GT\${TAG}.pc_2.png {1}_GT\${TAG}_Clusters.pc_2.png 
        # timeout 5h \${WIN} chromplot -f png {1}_GT\${TAG} {1}:500000-\${END} -m \${META} -g fog --numeric -p het -i 30 >> {1}_GT\${TAG}.log 2>&1
        # mv {1}_GT\${TAG}.het.png {1}_GT\${TAG}_Clusters.het.png 
    # fi

"
