#!/bin/bash

#SBATCH --job-name=ld_scan_outlier
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=esteiger@ucsc.edu
#SBATCH --output=logs/ld_scan_outlier_%A_%a.log

set -euo pipefail

LD_DIR=/private/groups/meyerlab/esteiger/variant_calling/ngsLD
WORKDIR=/private/groups/meyerlab/esteiger/selection/gemma/ld_outlier_windows
SCAN=/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/ld_scan_outlier_snps.py

cd "$WORKDIR"
mkdir -p linked_snp_positions logs

FILELIST=$1
R2_THRESHOLD=$2
SCAFFOLD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")

if [[ -z "$SCAFFOLD" ]]; then
    echo "No scaffold for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

OUT="linked_snp_positions/${SCAFFOLD}.linked.tsv"
if [[ -f "$OUT" ]]; then
    echo "Already done: ${SCAFFOLD}"
    exit 0
fi

LD_GZ="${LD_DIR}/${SCAFFOLD}_prune.ld.gz"
if [[ ! -f "$LD_GZ" ]]; then
    echo "LD file not found for ${SCAFFOLD}: ${LD_GZ}"
    exit 1
fi

python3 "$SCAN" \
    --scaffold "$SCAFFOLD" \
    --ld-gz "$LD_GZ" \
    --seeds outlier_snp_seeds.tsv \
    --r2-threshold "$R2_THRESHOLD" \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --out "${OUT}.tmp"
mv "${OUT}.tmp" "$OUT"

echo "Done: ${SCAFFOLD}"
