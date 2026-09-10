#!/bin/bash

#SBATCH --job-name=ngsLD
#SBATCH --partition=long    
#SBATCH --mail-type=BEGIN,END,FAIL          
#SBATCH --mail-user=esteiger@ucsc.edu 
#SBATCH --nodes=1    	
#SBATCH --cpus-per-task=120                 # 10 threads, 30 jobs
#SBATCH --mem=100gb                    # 5gb per job * 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs; medium is 
#SBATCH --output=serial_test_%j.log   # Standard output and error log

cd /private/groups/meyerlab/esteiger/variant_calling/ngsLD

find . -name "DEINVI_*.tsv.gz" | parallel -j 4 "./ngsLD.sh {} &> {}.out"
