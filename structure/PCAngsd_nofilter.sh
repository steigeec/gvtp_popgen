#!/bin/bash
# PCAngsd_nofilter.sh
# PCAngsd on all samples in a Beagle file.
# Usage: sbatch PCAngsd_nofilter.sh <beagle_prefix>
# Conda env: PCAngsd
#   conda create -n PCAngsd -c conda-forge python numpy scipy cython \
#     c-compiler cxx-compiler scikit-allel plotly python-kaleido
#   pip install -e <clone of https://github.com/Rosemeis/pcangsd>   # PCAngsd 1.36.0
#SBATCH --job-name=PCAngsd    
#SBATCH --output=PCAngsd.out  
#SBATCH --time=4:00:00            
#SBATCH --cpus-per-task=50                 
#SBATCH --nodes=1    	
#SBATCH --mem=50G                  
#SBATCH --partition=long       

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/structure/PCAngsd

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate PCAngsd || { echo "Conda environment activation failed"; exit 1; }
# Ensure system utilities are on PATH
export PATH="$PATH:/usr/bin"

SUBSET=$1
OUT_PATH="${BASE_DIR}/variant_calling/ngsLD"
BEAGLE=${OUT_PATH}/${SUBSET}.beagle.gz
PCANGSD="/private/home/esteiger/miniconda3/envs/PCAngsd/bin/pcangsd"

# MAF filter used because sites with MAF < 5% were already removed at the
# full-dataset scale, so the variant set is trusted; real variants found in only
# a small part of either subset could be interesting.
${PCANGSD} -b ${BEAGLE} -t 50 -o $SUBSET --maf 0.05 --iter 500 
