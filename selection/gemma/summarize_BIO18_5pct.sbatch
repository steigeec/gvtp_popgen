#!/bin/bash
# Summarize the BIO18 5% SNP-subset permutation null (100 permutations).

#SBATCH --job-name=bio18_5pct_summary
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=00:20:00
#SBATCH --output=/private/groups/meyerlab/esteiger/selection/gemma/logs/bio18_5pct_summary_%j.log

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/selection/gemma

Rscript ${BASE_DIR}/gvtp_popgen/selection/gemma/summarize_permutation_null.R \
  BIO18_5pct 100 BIO18
