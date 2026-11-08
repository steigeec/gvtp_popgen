#!/bin/bash
# GEMMA LMM permutation null for BIO18 on the 5% SNP subset (one permutation per array task).
# Conda env: gemma
#   conda create -n gemma -c bioconda gemma=0.98.3

#SBATCH --job-name=bio18_5pct_perm
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --array=1-100%25
#SBATCH --output=/private/groups/meyerlab/esteiger/selection/gemma/logs/bio18_5pct_perm_%A_%a.log

BASE_DIR=/private/groups/meyerlab/esteiger

eval "$(conda shell.bash hook)"
conda activate gemma

cd ${BASE_DIR}/selection/gemma

BFILE=postParalog_MAF_4x_GLimp_clean
KIN=output/kinship.cXX.txt
COV=covariates_intercept.txt
SNPS=snp_subset_5pct.txt
i=${SLURM_ARRAY_TASK_ID}

# perm_BIO18/perm_${i}.txt from make_permutations.R; permutations do not depend on the SNP subset.
gemma -bfile ${BFILE} -p perm_BIO18/perm_${i}.txt -k ${KIN} -lmm 2 -n 1 -c ${COV} -snps ${SNPS} -o BIO18_5pct_perm_${i}
