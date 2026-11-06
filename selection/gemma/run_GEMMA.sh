#!/bin/bash

#SBATCH --job-name=gemma_lmm
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=gemma_lmm_%A_%a.log
#SBATCH --array=1-16

eval "$(conda shell.bash hook)"
conda activate gemma

PCS="/private/groups/meyerlab/esteiger/selection/gemma/covariates_intercept.txt"
BFILE="/private/groups/meyerlab/esteiger/selection/gemma/postParalog_MAF_4x_GLimp_clean"
PHENO="/private/groups/meyerlab/esteiger/selection/gemma/env.txt"
KINSHIP="output/kinship.cXX.txt"
KINSHIP_DONE="output/kinship.done"

cd /private/groups/meyerlab/esteiger/selection/gemma

i=${SLURM_ARRAY_TASK_ID}

# Relatedness (kinship) matrix is shared across all 16 phenotypes, so it's
# computed once, gated to array task 1, instead of every task racing to
# write the same file. -n 1 only selects which individuals to include (any
# phenotype column with the same missingness pattern would do the same
# thing here) -- it does not restrict -gk to phenotype 1's values.
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

# run univariate LMM, kinship-corrected
# -k  relatedness matrix (shared across all phenotypes, computed above)
# -c  covariates, including an explicit intercept column (GEMMA does not
#     add one implicitly when -c is supplied -- see PC_env_creation.R)
# -lmm [num] option specifies which frequentist test to use, i.e. -lmm 1 performs
# Wald test, -lmm 2 performs likelihood ratio test, -lmm 3 performs score test, and -lmm 4
# performs all the three tests
# -n : which of the environmental variables do I want to use each time?
gemma -bfile ${BFILE} \
      -p ${PHENO} \
      -k ${KINSHIP} \
      -lmm 2 \
      -n $i \
      -c ${PCS} \
      -o pheno${i}_lmm_results
