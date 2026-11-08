#!/bin/bash

#SBATCH --job-name=FinishingGVTPmapping
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=50                    # number of CPUs
#SBATCH --mem=250gb                     # 5gb per job * 
#SBATCH --time=128:00:00               # long:168, medium: 12
#SBATCH --output=serial_test_%j.log   # Standard output and error log

# Map all samples in each trimmed-read batch directory with GVTP_mapping.sh.

BASE_DIR=/private/groups/meyerlab/esteiger

echo "Finishing gaviota tarplant mapping on Emerald"
date

subdirs=(
# batches already mapped are commented out
#"Guilliams_10089_240712A9"
#"Guilliams2_10089_240712B9"
#"Guilliams_10089_240723A9"
#"Guilliams_10089_240723B9"
"Guilliams_10089_240726A9"
#"K-M"
) 

for subdir in "${subdirs[@]}"; do
  cd ${BASE_DIR}/trimmed_reads/${subdir}
  outdir="${BASE_DIR}/mapping/${subdir}"
  find . -name "*_R1_clean.fastq.gz" | parallel -j 50 "${BASE_DIR}/gvtp_popgen/preprocessing/GVTP_mapping.sh {} $outdir 4 &> $outdir/{}.mapping.out"
done
