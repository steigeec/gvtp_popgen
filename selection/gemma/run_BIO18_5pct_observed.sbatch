#!/bin/bash
# GEMMA LMM for observed BIO18 on the 5% SNP subset used by the permutation array.
# Conda env: gemma
#   conda create -n gemma -c bioconda gemma=0.98.3

#SBATCH --job-name=bio18_5pct_obs
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/selection/gemma/logs/bio18_5pct_obs_%j.log

BASE_DIR=/private/groups/meyerlab/esteiger

eval "$(conda shell.bash hook)"
conda activate gemma

cd ${BASE_DIR}/selection/gemma

BFILE=postParalog_MAF_4x_GLimp_clean
KIN=output/kinship.cXX.txt
COV=covariates_intercept.txt
SNPS=snp_subset_5pct.txt

gemma -bfile ${BFILE} -p env.txt -k ${KIN} -lmm 2 -n 5 -c ${COV} -snps ${SNPS} -o BIO18_5pct_subset_observed
