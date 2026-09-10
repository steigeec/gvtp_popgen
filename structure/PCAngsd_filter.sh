#!/bin/bash


#SBATCH --job-name=PCAngsd    
#SBATCH --output=PCAngsd.out  
#SBATCH --time=4:00:00            
#SBATCH --cpus-per-task=50                 
#SBATCH --nodes=1    	
#SBATCH --mem=10G                  # CHECK AND REDUCE AS NEEDED
#SBATCH --partition=long       

cd /private/groups/meyerlab/esteiger/structure/PCAngsd

# Ensure Conda is initialized and activate environment
source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
# Activate PCAngsd environment
conda activate PCAngsd || { echo "Conda environment activation failed"; exit 1; }
# Explicitly add /usr/bin to PATH to ensure system utilities are available
export PATH="$PATH:/usr/bin"

BEAGLE_PREFIX=$1
OUT_NAME=${2:-$1}
OUT_PATH="/private/groups/meyerlab/esteiger/variant_calling/ngsLD"
BEAGLE=${OUT_PATH}/${BEAGLE_PREFIX}.beagle.gz
PCANGSD="/private/home/esteiger/miniconda3/envs/PCAngsd/bin/pcangsd"

# Run PCAngsd
${PCANGSD} -b ${BEAGLE} -t 50 -o ${OUT_NAME}_filter --maf 0.05 --iter 500 --filter ${OUT_NAME}_subset.txt


