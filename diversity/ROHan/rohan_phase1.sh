#!/bin/bash
# =============================================================================
# rohan_phase1.sh — SLURM array: ROHan main-run pass 1 (MCMC het estimation)
#
# One array task per remaining sample.  Task list read from rohan_tasks.tsv,
# rebuilt each time by launch_rohan.sh (only samples still missing a valid
# summary.txt are included).  Submit via launch_rohan.sh or manually:
#   sbatch --array=1-N%40 rohan_phase1.sh
#
# Skip check validates the existing .hEst.gz with `gzip -t`, not just
# non-empty (-s) — a file left behind by a job that was killed mid-write
# (e.g. hit the SLURM time limit) is non-empty but truncated, and silently
# treating it as "done" would feed pass 2 a corrupt input.
# =============================================================================

#SBATCH --job-name=rohan_p1
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --time=3-00:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/logs/rohan_phase1_%A_%a.log

set -euo pipefail

ROHAN=/private/groups/meyerlab/esteiger/software/rohan/bin/rohan
REF=/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked
SCAFFOLDS=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/scaffolds_1mb_names.txt
TASK_FILE=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/rohan_tasks.tsv
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
