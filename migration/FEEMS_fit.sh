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
# Conda env: feems_e
#   conda create -n feems_e -c conda-forge -c bioconda python=3.8.10 \
#     numpy=1.22.3 scipy=1.5.0 scikit-learn=0.23.1 matplotlib=3.2.2 \
#     networkx=2.4 pyproj=2.6.1 shapely=1.7.1 fiona cartopy=0.18.0 \
#     suitesparse=5.7.2 scikit-sparse=0.4.4 pandas-plink msprime=1.0.0 \
#     statsmodels=0.14.0
#   pip install git+https://github.com/NovembreLab/feems@c08a179   # FEEMS 1.0.1

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/migration/FEEMS

N=${SLURM_ARRAY_TASK_ID}

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate feems_e || { echo "Failed to activate Conda environment"; exit 1; }
echo "conda path is $CONDA_PREFIX"

$CONDA_PREFIX/bin/python ${BASE_DIR}/gvtp_popgen/migration/FEEMS_fit.py year${N}_strauss_maf
