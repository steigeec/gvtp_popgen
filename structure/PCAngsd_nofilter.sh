#!/bin/bash

#SBATCH --job-name=PCAngsd    
#SBATCH --output=PCAngsd.out  
#SBATCH --time=4:00:00            
#SBATCH --cpus-per-task=50                 
#SBATCH --nodes=1    	
#SBATCH --mem=50G                  
#SBATCH --partition=long       

cd /private/groups/meyerlab/esteiger/structure/PCAngsd

# Ensure Conda is initialized and activate environment
source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
# Activate PCAngsd environment
conda activate PCAngsd || { echo "Conda environment activation failed"; exit 1; }
# Explicitly add /usr/bin to PATH to ensure system utilities are available
export PATH="$PATH:/usr/bin"

SUBSET=$1
OUT_PATH="/private/groups/meyerlab/esteiger/variant_calling/ngsLD"
BEAGLE=${OUT_PATH}/${SUBSET}.beagle.gz
PCANGSD="/private/home/esteiger/miniconda3/envs/PCAngsd/bin/pcangsd"

# Use the maf filter because we already threw out sites with a MAF<5% at the full dataset scale
#       so we feel confident in the variant set that we have
#       real variants only occuring in a small portion of the data for either of these subsets could actually be interesting

# Run PCAngsd
${PCANGSD} -b ${BEAGLE} -t 50 -o $SUBSET --maf 0.05 --iter 500 
