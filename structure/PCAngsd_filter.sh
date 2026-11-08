#!/bin/bash
# PCAngsd_filter.sh
# PCAngsd on a subset of samples selected by <out_name>_subset.txt (from run_PCAngsd.sh).
# Usage: sbatch PCAngsd_filter.sh <beagle_prefix> [out_name]
# Conda env: PCAngsd
#   conda create -n PCAngsd -c conda-forge python numpy scipy cython \
#     c-compiler cxx-compiler scikit-allel plotly python-kaleido
#   pip install -e <clone of https://github.com/Rosemeis/pcangsd>   # PCAngsd 1.36.0
#SBATCH --job-name=PCAngsd    
#SBATCH --output=PCAngsd.out  
#SBATCH --time=4:00:00            
#SBATCH --cpus-per-task=50                 
#SBATCH --nodes=1    	
#SBATCH --mem=10G                  # CHECK AND REDUCE AS NEEDED
#SBATCH --partition=long       

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/structure/PCAngsd

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate PCAngsd || { echo "Conda environment activation failed"; exit 1; }
# Ensure system utilities are on PATH
export PATH="$PATH:/usr/bin"

BEAGLE_PREFIX=$1
OUT_NAME=${2:-$1}
OUT_PATH="${BASE_DIR}/variant_calling/ngsLD"
BEAGLE=${OUT_PATH}/${BEAGLE_PREFIX}.beagle.gz
PCANGSD="/private/home/esteiger/miniconda3/envs/PCAngsd/bin/pcangsd"

${PCANGSD} -b ${BEAGLE} -t 50 -o ${OUT_NAME}_filter --maf 0.05 --iter 500 --filter ${OUT_NAME}_subset.txt

