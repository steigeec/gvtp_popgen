#!/bin/bash

#SBATCH --job-name=bio18_5pct_obs
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/selection/gemma/logs/bio18_5pct_obs_%j.log

eval "$(conda shell.bash hook)"
conda activate gemma

cd /private/groups/meyerlab/esteiger/selection/gemma

BFILE=postParalog_MAF_4x_GLimp_clean
KIN=output/kinship.cXX.txt
COV=covariates_intercept.txt
SNPS=snp_subset_5pct.txt

# Observed (unpermuted) BIO18 result on the same 5% SNP subset used for the
# permutation array, so observed vs. null is an apples-to-apples comparison.
gemma -bfile ${BFILE} -p env.txt -k ${KIN} -lmm 2 -n 5 -c ${COV} -snps ${SNPS} -o BIO18_5pct_subset_observed
