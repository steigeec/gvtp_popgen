#!/bin/bash --login

#SBATCH --job-name=FRAME_lamb
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --array=1-3
#SBATCH --output=/private/groups/meyerlab/esteiger/migration/FEEMS/FRAME_lambda_%a.log

# Refit the pooled FRAME model (strauss_all_maf) at fixed, more strongly
# smoothed lambda_m values, one per array task, for frame_wind_asymmetry.R.
# Conda env: frame_e
#   conda create -n frame_e -c conda-forge python=3.11.9 numpy=1.26.4 \
#     scipy=1.11.4 scikit-learn=1.5.1 matplotlib=3.10 networkx=3.4.2 \
#     pyproj=3.6.1 shapely=2.0.6 fiona=1.10.1 cartopy=0.24.1 \
#     "blas=*=mkl" mkl-service
#   pip install git+https://github.com/NovembreLab/frame@98abf4b   # FRAME 1.0.0
#   (run with PYTHONNOUSERSITE=1 so packages in ~/.local are not picked up)

BASE_DIR=/private/groups/meyerlab/esteiger

LAMBDAS=(0.0379 0.336 2.98)

cd ${BASE_DIR}/migration/FEEMS
source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate frame_e || { echo "Failed to activate Conda environment"; exit 1; }
export PYTHONNOUSERSITE=1

"$CONDA_PREFIX/bin/python" ${BASE_DIR}/gvtp_popgen/migration/FRAME_fit.py \
    strauss_all_maf "${LAMBDAS[$((SLURM_ARRAY_TASK_ID - 1))]}"
