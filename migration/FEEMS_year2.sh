#!/bin/bash --login

#SBATCH --job-name=FEEMS_y2
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=300gb
#SBATCH --time=168:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/migration/FEEMS/FEEMS_year2.log

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/migration/FEEMS

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate feems_e || { echo "Failed to activate Conda environment"; exit 1; }
echo "conda path is $CONDA_PREFIX"

$CONDA_PREFIX/bin/python ${BASE_DIR}/gvtp_popgen/migration/FEEMS.py year2_strauss_maf
