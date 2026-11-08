#!/bin/bash

#SBATCH --job-name=Merging_BAMs
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=30                  # number of CPUs
#SBATCH --mem=210gb                     # 7gb per job * 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=serial_test_%j.log   # Standard output and error log

# Run merge_bams.sh for every sample ID in sampleIDs.txt.

export BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/mapping
input="${BASE_DIR}/sampleIDs.txt"

ulimit -m 40000

parallel -a $input -j 15 '${BASE_DIR}/gvtp_popgen/preprocessing/merge_bams.sh {} > ${BASE_DIR}/merged/{}_merge.log 2>&1'
