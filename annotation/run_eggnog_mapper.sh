#!/bin/bash

#SBATCH --job-name=eggnog_mapper
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=48:00:00
#SBATCH --output=eggnog_mapper_%j.log

# eggNOG-mapper on the reference CDS (orthology-based GO annotation); reshape
# the output with convert_eggnog_to_geneid_go.py.
# Conda env: eggnog-mapper
#   built by annotation/setup_eggnog_mapper.sh (plus its database download):
#   conda create -n eggnog-mapper -c bioconda -c conda-forge eggnog-mapper=2.1.15
# Usage: sbatch run_eggnog_mapper.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

CONDA_BASE=/private/home/esteiger/miniconda3
ENV_NAME=eggnog-mapper
DATA_DIR=${BASE_DIR}/software/eggnog-mapper-data

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

CDS=${BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.cds
OUTDIR=${BASE_DIR}/genome/annotation
OUT_PREFIX=gvtp_eggnog_Di2_viridiplantae

mkdir -p "${OUTDIR}/tmp"
cd "${OUTDIR}"

# --itype CDS: input is nucleotide; emapper translates before DIAMOND
# --tax_scope Viridiplantae: restrict GO transfer to green-plant orthologs
emapper.py \
    -i "${CDS}" \
    --itype CDS \
    -m diamond \
    --sensmode sensitive \
    --tax_scope Viridiplantae \
    --go_evidence non-electronic \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    --data_dir "${DATA_DIR}" \
    --temp_dir "${OUTDIR}/tmp" \
    --output "${OUT_PREFIX}" \
    --output_dir "${OUTDIR}" \
    --override

echo "Done. Main output: ${OUTDIR}/${OUT_PREFIX}.emapper.annotations"
echo "Next: python3 convert_eggnog_to_geneid_go.py ${OUTDIR}/${OUT_PREFIX}.emapper.annotations"
