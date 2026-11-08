#!/bin/bash

#SBATCH --job-name=launch_heatmapLD
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=1                 # 10 threads, 30 jobs
#SBATCH --mem=200gb                     # 5gb per job * 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs; medium is 
#SBATCH --output=serial_test_%j.log   # Standard output and error log

# Run heatmap_processed.sh on the DEINVI_00007 LD files.

export BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/variant_calling/ngsLD

find . -name 'DEINVI_00007*_heatmap_processed.ld' | parallel -j 1 '
	${BASE_DIR}/gvtp_popgen/variant_calling_filtering/heatmap_processed.sh $(basename {})
'

