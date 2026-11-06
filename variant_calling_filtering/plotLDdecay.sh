#!/bin/bash


#SBATCH --job-name=plotLDdecay
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=1                  # 8 threads per job, 23 jobs
#SBATCH --mem=1000gb                     # 5gb per job * 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=serial_test_%j.log   # Standard output and error log

cd /private/groups/meyerlab/esteiger/variant_calling/ngsLD
SCR="/private/groups/meyerlab/esteiger/software/ngsLD/scripts/fit_LDdecay.R" 
R_LIBS="/private/groups/meyerlab/esteiger/software/R_libs" Rscript --vanilla --slave $SCR --ld_files myLDfiles.txt --out angsd_LD.pdf --fit_level 10 --out angsd_LD.pdf --max_kb_dist 300 --min_maf 0.05 --n_ind 842


