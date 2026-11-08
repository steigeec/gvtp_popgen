#!/bin/bash
# GEMMA univariate LMM (likelihood ratio test) for each of the 16 environmental variables in env.txt, one per array task.
# Conda env: gemma
#   conda create -n gemma -c bioconda gemma=0.98.3

#SBATCH --job-name=gemma_lmm
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=gemma_lmm_%A_%a.log
#SBATCH --array=1-16

BASE_DIR=/private/groups/meyerlab/esteiger

eval "$(conda shell.bash hook)"
conda activate gemma

PCS="${BASE_DIR}/selection/gemma/covariates_intercept.txt"
BFILE="${BASE_DIR}/selection/gemma/postParalog_MAF_4x_GLimp_clean"
PHENO="${BASE_DIR}/selection/gemma/env.txt"
KINSHIP="output/kinship.cXX.txt"
KINSHIP_DONE="output/kinship.done"

cd ${BASE_DIR}/selection/gemma

i=${SLURM_ARRAY_TASK_ID}

# Kinship matrix is computed once by task 1; -n 1 only sets which individuals are included.
if [ "$i" -eq 1 ]; then
    if [ -f "$KINSHIP_DONE" ]; then
        echo "Kinship matrix already computed, skipping."
    else
        echo "Computing kinship matrix..."
        gemma -bfile ${BFILE} -p ${PHENO} -n 1 -gk 1 -o kinship
        touch "$KINSHIP_DONE"
    fi
else
    echo "Waiting for kinship matrix from array task 1..."
    waited=0
    while [ ! -f "$KINSHIP_DONE" ]; do
        if [ "$waited" -ge 5400 ]; then
            echo "Timed out after 90 min waiting for ${KINSHIP_DONE} -- task 1 may have failed."
            exit 1
        fi
        sleep 30
        waited=$((waited + 30))
    done
    echo "Kinship matrix ready."
fi

# -c covariates include an explicit intercept column (GEMMA does not add one when -c is given).
# -lmm 2 = likelihood ratio test; -n = column of env.txt.
gemma -bfile ${BFILE} \
      -p ${PHENO} \
      -k ${KINSHIP} \
      -lmm 2 \
      -n $i \
      -c ${PCS} \
      -o pheno${i}_lmm_results
