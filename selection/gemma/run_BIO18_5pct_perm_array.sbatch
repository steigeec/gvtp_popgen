#!/bin/bash

#SBATCH --job-name=bio18_5pct_perm
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --array=1-100%25
#SBATCH --output=/private/groups/meyerlab/esteiger/selection/gemma/logs/bio18_5pct_perm_%A_%a.log

eval "$(conda shell.bash hook)"
conda activate gemma

cd /private/groups/meyerlab/esteiger/selection/gemma

BFILE=postParalog_MAF_4x_GLimp_clean
KIN=output/kinship.cXX.txt
COV=covariates_intercept.txt
SNPS=snp_subset_5pct.txt
i=${SLURM_ARRAY_TASK_ID}

# perm_BIO18/perm_${i}.txt was already generated (make_permutations.R) for the
# earlier 1% subset run — the permutation itself only shuffles phenotype-
# individual pairing, so it's independent of which SNP subset is analyzed and
# is safe to reuse here unchanged.
gemma -bfile ${BFILE} -p perm_BIO18/perm_${i}.txt -k ${KIN} -lmm 2 -n 1 -c ${COV} -snps ${SNPS} -o BIO18_5pct_perm_${i}
