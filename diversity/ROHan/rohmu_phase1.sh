#!/bin/bash
# rohmu_phase1.sh
# SLURM array, ROHan pass 1 (MCMC) for rohmu calibration, one sample at one
# depth per task from rohmu_calibration/phase1_tasks.tsv. Submit via
# estimate_rohmu.sh or:
#   sbatch --array=1-N rohmu_phase1.sh

#SBATCH --job-name=rohmu_p1
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=168:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/rohmu_calibration/logs/phase1_%A_%a.log

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

ROHAN=${BASE_DIR}/software/rohan/bin/rohan
REF=${BASE_DIR}/genome/DEINVI_genome.fasta.masked
SCAFFOLDS=${BASE_DIR}/diversity/ROH/ROHan/scaffolds_1mb_names.txt
TASK_FILE=${BASE_DIR}/diversity/ROH/ROHan/rohmu_calibration/phase1_tasks.tsv

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_FILE")
IFS=$'\t' read -r prefix bam chains <<< "$LINE"

echo "Task ${SLURM_ARRAY_TASK_ID}: prefix=$prefix bam=$bam chains=$chains"

if [[ -s "${prefix}.hEst.gz" ]]; then
    echo "Pass 1 already complete, skipping."
    exit 0
fi

"$ROHAN" \
    --tstv 1.86 --size 100000 --chains "$chains" --nohmm \
    --auto "$SCAFFOLDS" \
    -o "$prefix" "$REF" "$bam" \
    2>> "${prefix}.log"

echo "Pass 1 complete."
