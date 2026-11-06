#!/bin/bash

#SBATCH --job-name=delly
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=200                  
#SBATCH --mem=500gb                     # 1gb per task with winPCA 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=DELLY_%j.log   # Standard output and error log


cd /private/groups/meyerlab/esteiger/inversions

# Delly primarily parallelizes on the sample level. 
# Hence, OMP_NUM_THREADS should be always <= to the number of input samples.
# The number of samples is 293
export OMP_NUM_THREADS=200

DELLY="/private/groups/meyerlab/esteiger/software/delly/src/delly"
OUT="structural_variants.bcf"
GENOME="/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked"

BAMLIST="/private/groups/meyerlab/esteiger/mapping/merged/bamlist_4x.txt"
BAM=$(sed 's|^|/private/groups/meyerlab/esteiger/mapping/merged/|' "$BAMLIST" | tr '\n' ' ' | sed 's/ *$//')


${DELLY} call -o ${OUT} -g ${GENOME} ${BAM}

