#!/bin/bash

#SBATCH --job-name=winPCA
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=200    # also use --cpus-per-task=100 when actually calculating tsvs?             
#SBATCH --mem=200gb                     # 1gb per task with winPCA 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=serial_test_%j.log   # Standard output and error log

cd /private/groups/meyerlab/esteiger/inversions/all_localities

FAI='/private/groups/meyerlab/esteiger/inversions/filtered_scaffolds.fai'
BEAGLE='/private/groups/meyerlab/esteiger/variant_calling/angsd/paralog_filtered/four_coverage_beforeMAF.beagle'
VCF='/private/groups/meyerlab/esteiger/variant_calling/angsd/paralog_filtered/four_coverage_beforeMAF.vcf'
SUBSET='/private/groups/meyerlab/esteiger/mapping/merged/bamlist_8x.txt'
WIN='/private/groups/meyerlab/esteiger/software/winpca/winpca'
CONDA_SH=/private/groups/meyerlab/esteiger/software/miniconda3/etc/profile.d/conda.sh
META='/private/groups/meyerlab/esteiger/inversions/metadata.tsv'

COLORS='Strauss:009e73,CDFW:e69f00,Dangermond:56b4e9,VAFB_Arguello:d55e00,VAFB_Lion_Head:f0e442,Point_Conception:0072b2'

export FAI
export COLORS
export CONDA_SH
export WIN
export BEAGLE
export VCF
export SUBSET
export META

# Before the parallel command, set up the full PATH
export PATH=/usr/bin:/bin:/private/groups/meyerlab/esteiger/software/miniconda3/bin:$PATH


cat ${FAI} | parallel -j 5  --colsep '\t' --env META --env COLORS --env PATH --env CONDA_SH --env WIN --env BEAGLE "

    WINDOW=1000000
    INCREMENT=10000
    HALF=\$((WINDOW/2))
    END=\$(({2} - HALF))
    # Source conda setup in each subshell
    source \${CONDA_SH}
    # Create log file for a given chrom
    LOG={1}.log    

     echo \"DEBUG: Processing scaffold \$SCAFFOLD with length \$LENGTH\" >> \${SCAFFOLD}.log 2>&1

    # Skip iteration if column 2 is not > 800,000
    if [ {2} -le 900000 ]; then
        echo 'Skipping {1} because {2} < 900,000' >> {1}.log 2>&1
        exit 0
    fi

    # Skip if there's already an out file
    if [ -f \${LOG} ]; then
        echo 'skipping because log file already exists' >> {1}.log 2>&1
        exit 0
    fi

    # GL-based PCA, standard inrement of 10k
    if [ ! -f {1}_GL.stat.tsv.gz ]; then
        \${WIN} pca -v GL -t 40 -w 1000000 --increment 10000 \${BEAGLE} {1}:1-{2} {1}_GL >> {1}.log 2>&1
    fi

    # Ensure the stat file exists before using it
    if [ ! -f {1}_GL.stat.tsv.gz ]; then
        echo 'Error: {1}_GL.stat.tsv.gz not found!' >> {1}.log 2>&1
        exit 1
    fi

    # PC plots
    if [ ! -f {1}_GL.pc_1.png ] && [ -f {1}_GL.stat.tsv.gz ]; then
        # When I do use interval, I have to call like this        \${INTERVAL}
        # -c \${COLORS}
        #timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_location.pc_1.png
        timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m \${META} -g fog --numeric >> {1}.log 2>&1
	mv {1}_GL.pc_1.png {1}_GL_fog.pc_1.png
        #timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m metadata.tsv -g deficit --numeric >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_deficit.pc_1.png
    fi

    if [ ! -f {1}_GL.pc_1.png ] && [ -f {1}_GL.stat.tsv.gz ]; then
        echo 'attempting for interval of 5' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GL_location {1}:500000-\${END} -m metadata.tsv -g Clusters -i 5 -c \${COLORS} >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_location.pc_1.png
        timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m \${META} -g fog -i 5 --numeric >> {1}.log 2>&1
	mv {1}_GL.pc_1.png {1}_GL_fog.pc_1.png
        #timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m metadata.tsv -g deficit -i 5 --numeric >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_deficit.pc_1.png
    fi

    if [ ! -f {1}_GL.pc_1.png ] && [ -f {1}_GL.stat.tsv.gz ]; then
        echo 'attempting for interval of 10' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GL_location {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -i 10 >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_location.pc_1.png
        timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m \${META} -g fog -i 10 --numeric >> {1}.log 2>&1
	mv {1}_GL.pc_1.png {1}_GL_fog.pc_1.png
        #timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m metadata.tsv -g deficit -i 10 --numeric >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_deficit.pc_1.png
    fi

    if [ ! -f {1}_GL.pc_1.png ] && [ -f {1}_GL.stat.tsv.gz ]; then
        echo 'attempting for interval of 15' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GL_location {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -i 15 >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_location.pc_1.png
        timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m \${META} -g fog -i 15 --numeric >> {1}.log 2>&1
	mv {1}_GL.pc_1.png {1}_GL_fog.pc_1.png
        #timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m metadata.tsv -g deficit -i 15 --numeric >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_deficit.pc_1.png
    fi

    if [ ! -f {1}_GL.pc_1.png ] && [ -f {1}_GL.stat.tsv.gz ]; then
        echo 'attempting for interval of 20' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GL_location {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -i 20 >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_location.pc_1.png
        timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m \${META} -g fog -i 20 --numeric >> {1}.log 2>&1
	mv {1}_GL.pc_1.png {1}_GL_fog.pc_1.png
        #timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m metadata.tsv -g deficit -i 20 --numeric >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_deficit.pc_1.png
    fi

    if [ ! -f {1}_GL.pc_1.png ] && [ -f {1}_GL.stat.tsv.gz ]; then
        echo 'attempting for interval of 25' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GL_location {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -i 25 >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_location.pc_1.png
        timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m \${META} -g fog -i 25 --numeric >> {1}.log 2>&1
	mv {1}_GL.pc_1.png {1}_GL_fog.pc_1.png        
	#timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m metadata.tsv -g deficit -i 25 --numeric >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_deficit.pc_1.png
    fi

    if [ ! -f {1}_GL.pc_1.png ] && [ -f {1}_GL.stat.tsv.gz ]; then
        echo 'attempting for interval of 30' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GL_location {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -i 30 >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_location.pc_1.png
        timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m \${META} -g fog -i 30 --numeric >> {1}.log 2>&1
	mv {1}_GL.pc_1.png {1}_GL_fog.pc_1.png
        #timeout 5h \${WIN} chromplot -f png {1}_GL {1}:500000-\${END} -m metadata.tsv -g deficit -i 30 --numeric >> {1}.log 2>&1
	#mv {1}_GL.pc_1.png {1}_GL_deficit.pc_1.png
    fi
    
    exit 1
     
    #########################################

    # GT-based PCA
    
    # Restrict to 8x+ only... -s \${SUBSET}
    if [ ! -f {1}_GT.stat.tsv.gz ]; then
        \${WIN} pca -v GT -t 40 -w 1000000 --increment 10000 \${VCF} {1}:1-{2} {1}_GT >> {1}.log 2>&1
    fi

    # Ensure the stat file exists before using it
    if [ ! -f {1}_GT.stat.tsv.gz ]; then
        echo 'Error: {1}_GT.stat.tsv.gz for genotype-based work not found!' >> {1}.log 2>&1
        exit 1
    fi

    # GRADUALLY TRY LARGER AND LARGER INTERVALS TO ACHIEVE SUCCESSFUL PLOTTING... 

    # PC plots
    if [ ! -f {1}_GT_Clusters.pc_1.png ] && [ -f {1}_GT.stat.tsv.gz ]; then
        echo 'attempting for no interval' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric >> {1}.log 2>&1
        #mv {1}_GT.pc_1.png {1}_GT_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -p PC2 >> {1}.log 2>&1
        #mv {1}_GT.pc_2.png {1}_GT_Clusters.pc_2.png 
        timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -p het >> {1}.log 2>&1
        mv {1}_GT.het.png {1}_GT_Clusters.het.png 
    fi

    if [ ! -f {1}_GT_Clusters.pc_1.png ] && [ -f {1}_GT.stat.tsv.gz ]; then
        echo 'attempting for interval of 5' 
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -i 5 >> {1}.log 2>&1
        #mv {1}_GT.pc_1.png {1}_GT_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -p PC2 -i 5>> {1}.log 2>&1
        #mv {1}_GT.pc_2.png {1}_GT_Clusters.pc_2.png 
        timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -p het -i 5 >> {1}.log 2>&1
        mv {1}_GT.het.png {1}_GT_Clusters.het.png 
    fi

    if [ ! -f {1}_GT_Clusters.pc_1.png ] && [ -f {1}_GT.stat.tsv.gz ]; then
        echo 'attempting for interval of 10' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -i 10 >> {1}.log 2>&1
        #mv {1}_GT.pc_1.png {1}_GT_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -p PC2 -i 10 >> {1}.log 2>&1
        #mv {1}_GT.pc_2.png {1}_GT_Clusters.pc_2.png 
        timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -p het -i 10 >> {1}.log 2>&1
        mv {1}_GT.het.png {1}_GT_Clusters.het.png 
    fi

    if [ ! -f {1}_GT_Clusters.pc_1.png ] && [ -f {1}_GT.stat.tsv.gz ]; then
        echo 'attempting for interval of 15' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -i 15 >> {1}.log 2>&1
        #mv {1}_GT.pc_1.png {1}_GT_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -p PC2 -i 15>> {1}.log 2>&1
        #mv {1}_GT.pc_2.png {1}_GT_Clusters.pc_2.png 
        timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -p het -i 15 >> {1}.log 2>&1
        mv {1}_GT.het.png {1}_GT_Clusters.het.png 
    fi

    if [ ! -f {1}_GT_Clusters.pc_1.png ] && [ -f {1}_GT.stat.tsv.gz ]; then
        echo 'attempting for interval of 20' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -i 20 >> {1}.log 2>&1
        #mv {1}_GT.pc_1.png {1}_GT_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -p PC2 -i 20 >> {1}.log 2>&1
        #mv {1}_GT.pc_2.png {1}_GT_Clusters.pc_2.png 
        timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -p het -i 20 >> {1}.log 2>&1
        mv {1}_GT.het.png {1}_GT_Clusters.het.png 
    fi

    if [ ! -f {1}_GT_Clusters.pc_1.png ] && [ -f {1}_GT.stat.tsv.gz ]; then
        echo 'attempting for interval of 25' >> {1}.log 2>&1
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -i 25 >> {1}.log 2>&1
        #mv {1}_GT.pc_1.png {1}_GT_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -p PC2 -i 25>> {1}.log 2>&1
        #mv {1}_GT.pc_2.png {1}_GT_Clusters.pc_2.png 
        timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -p het -i 25 >> {1}.log 2>&1
        mv {1}_GT.het.png {1}_GT_Clusters.het.png 
    fi

    if [ ! -f {1}_GT_Clusters.pc_1.png ] && [ -f {1}_GT.stat.tsv.gz ]; then
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -i 30 >> {1}.log 2>&1
        #mv {1}_GT.pc_1.png {1}_GT_Clusters.pc_1.png 
        #timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g Clusters -c \${COLORS} -p PC2 -i 30>> {1}.log 2>&1
        #mv {1}_GT.pc_2.png {1}_GT_Clusters.pc_2.png 
        timeout 5h \${WIN} chromplot -f png {1}_GT {1}:500000-\${END} -m metadata.tsv -g fog --numeric -p het -i 30 >> {1}.log 2>&1
        mv {1}_GT.het.png {1}_GT_Clusters.het.png 
    fi

"
