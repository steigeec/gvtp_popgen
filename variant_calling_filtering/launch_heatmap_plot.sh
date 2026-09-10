#!/bin/bash

#SBATCH --job-name=launch_heatmapLD
#SBATCH --partition=long    
#SBATCH --mail-type=BEGIN,END,FAIL          
#SBATCH --mail-user=esteiger@ucsc.edu 
#SBATCH --nodes=1    	
#SBATCH --ntasks=1                 # 10 threads, 30 jobs
#SBATCH --mem=200gb                     # 5gb per job * 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs; medium is 
#SBATCH --output=serial_test_%j.log   # Standard output and error log

cd /private/groups/meyerlab/esteiger/variant_calling/ngsLD

find . -name 'DEINVI_00007*_heatmap_processed.ld' | parallel -j 1 '
	./heatmap_processed.R $(basename {})
'

