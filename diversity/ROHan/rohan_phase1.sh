#!/bin/bash
# rohan_phase1.sh
# SLURM array, ROHan pass 1 (MCMC heterozygosity estimation, --nohmm), one task
# per sample in rohan_tasks.tsv. Submit via launch_rohan.sh or:
#   sbatch --array=1-N%40 rohan_phase1.sh
#
# Existing hEst.gz files are checked with gzip -t, since a job killed mid-write
# leaves a non-empty but truncated file.

#SBATCH --job-name=rohan_p1
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --time=3-00:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/logs/rohan_phase1_%A_%a.log

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

ROHAN=${BASE_DIR}/software/rohan/bin/rohan
REF=${BASE_DIR}/genome/DEINVI_genome.fasta.masked
SCAFFOLDS=${BASE_DIR}/diversity/ROH/ROHan/scaffolds_1mb_names.txt
TASK_FILE=${BASE_DIR}/diversity/ROH/ROHan/rohan_tasks.tsv
CHAINS=10000
THREADS=8

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_FILE")
IFS=$'\t' read -r prefix bam <<< "$LINE"

echo "Task ${SLURM_ARRAY_TASK_ID}: prefix=$prefix bam=$bam"

if [[ -s "${prefix}.hEst.gz" ]] && gzip -t "${prefix}.hEst.gz" 2>/dev/null; then
    echo "Pass 1 already complete and valid, skipping."
    exit 0
fi

"$ROHAN" \
    --tstv 1.86 --size 100000 --chains "$CHAINS" --nohmm \
    -t "$THREADS" \
    --auto "$SCAFFOLDS" \
    -o "$prefix" "$REF" "$bam" \
    2>> "${prefix}.log"

echo "Pass 1 complete."
