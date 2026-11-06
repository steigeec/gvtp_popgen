#!/bin/bash --login

#SBATCH --job-name=FEEMS_fit
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=150gb
#SBATCH --time=168:00:00
#SBATCH --array=1-3
#SBATCH --output=/private/groups/meyerlab/esteiger/migration/FEEMS/FEEMS_cv_year%a.log

# FEEMS per sampling year (year1 = 2020, year2 = 2021, year3 = 2022).
# Lambda is cross-validated; append a value as a second argument to fix it.
# Inputs (year{N}_strauss_maf.*, grid.csv, edges.csv, outer.csv) come from FEEMS.sh.

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/migration/FEEMS

N=${SLURM_ARRAY_TASK_ID}

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate feems_e || { echo "Failed to activate Conda environment"; exit 1; }
echo "conda path is $CONDA_PREFIX"

$CONDA_PREFIX/bin/python ${BASE_DIR}/gvtp_popgen/migration/FEEMS_fit.py year${N}_strauss_maf
