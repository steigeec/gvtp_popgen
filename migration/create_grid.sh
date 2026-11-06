#!/bin/bash

#SBATCH --job-name=create_grid
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=1                  
#SBATCH --mem=2000gb                     # 5gb per job * 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=/private/groups/meyerlab/esteiger/migration/FEEMS/create_grid_%j.log   # Standard output and error log

BASE_DIR=/private/groups/meyerlab/esteiger

export R_LIBS_USER="/public/groups/meyerlab/esteiger/R_libs"

cd ${BASE_DIR}/migration/FEEMS

Rscript ${BASE_DIR}/gvtp_popgen/migration/create_grid.R
