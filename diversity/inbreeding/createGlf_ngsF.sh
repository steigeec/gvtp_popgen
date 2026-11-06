#!/bin/bash

#SBATCH --job-name=createGlf_ngsF
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=9                  # 8 threads per job, 23 jobs
#SBATCH --mem=150gb                    # bumped from 50gb: chunks 10-14 OOM'd at 50gb (MaxRSS ~52gb)
#SBATCH --time=48:00:00               # bumped down from 168h now that -rf regions file is deduplicated (~5h/chunk observed)
#SBATCH --output=serial_test_%j.log   # Standard output and error log


GENOME="/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked"
SITES="/private/groups/meyerlab/esteiger/variant_calling/angsd/subsample.sites"
chroms="/private/groups/meyerlab/esteiger/variant_calling/angsd/subsample.chrom"
BAMS=$1
OUT=/private/groups/meyerlab/esteiger/diversity/ngsF/${BAMS}

exec >${OUT}.log 2>&1

cd /private/groups/meyerlab/esteiger/mapping/merged

angsd -bam ${OUT} -sites $SITES -rf $chroms -ref $GENOME -P 9 -out $OUT \
-GL 2 -minMapQ 20 -minQ 20 -doGlf 3 -doMaf 1 -doCounts 1 -doHWE 1 -doMajorMinor 4 \
-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 0 -skipTriallelic 1  
