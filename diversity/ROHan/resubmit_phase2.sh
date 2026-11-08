#!/bin/bash
# resubmit_phase2.sh
# Resubmit a phase 2 array over all samples to pick up those whose phase 1 has
# finished. Safe to run while phase 1 is running, since phase 2 only reads
# hEst.gz and each task skips samples that are not ready or already done.
#
# Uses rohan_tasks_all.tsv (all 423 samples), kept separate from
# rohan_tasks.tsv because phase 1 reads that file by line number.
#
# Usage: bash resubmit_phase2.sh   (run on a schedule, e.g. cron)

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

OUTDIR=${BASE_DIR}/diversity/ROH/ROHan
SCRIPT_DIR=${BASE_DIR}/gvtp_popgen/diversity/ROHan
TASK_FILE=$OUTDIR/rohan_tasks_all.tsv
MAX_CONCURRENT_P2=40

[[ -s "$TASK_FILE" ]] || {
    echo "$(date): ERROR: $TASK_FILE not found or empty — nothing to submit." >&2
    exit 1
}

n=$(wc -l < "$TASK_FILE")
if (( n == 0 )); then
    echo "$(date): rohan_tasks_all.tsv is empty — nothing to submit."
    exit 0
fi

echo "$(date): Resubmitting phase2 array: $n tasks, ${MAX_CONCURRENT_P2} concurrent max"
P2=$(sbatch --parsable --array=1-${n}%${MAX_CONCURRENT_P2} "$SCRIPT_DIR/rohan_phase2.sh")
echo "$(date): Phase 2 job ID: $P2"
