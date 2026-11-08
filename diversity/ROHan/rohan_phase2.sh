#!/bin/bash
# rohan_phase2.sh
# SLURM array, ROHan pass 2 (HMM ROH calling), one task per sample in
# rohan_tasks_all.tsv (all 423 samples). rohmu is read from chosen_rohmu.txt.
#
# Uses its own task file because rohan_phase1.sh reads rohan_tasks.tsv by line
# number, so that file cannot change while phase 1 tasks are pending.
# Submitted without a dependency on phase 1; tasks whose hEst.gz is missing or
# invalid exit 0 and are picked up by resubmit_phase2.sh.
#
# Resources: across 466 HMM runs at rohmu=1e-5, max elapsed was ~17 min and max
# RSS ~10 MB. At higher rohmu a minority of samples take longer, hence the 2 h
# limit on the medium partition.

#SBATCH --job-name=rohan_p2
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=512mb
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/logs/rohan_phase2_%A_%a.log

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

ROHAN=${BASE_DIR}/software/rohan/bin/rohan
REF=${BASE_DIR}/genome/DEINVI_genome.fasta.masked
OUTDIR=${BASE_DIR}/diversity/ROH/ROHan
SCAFFOLDS=$OUTDIR/scaffolds_1mb_names.txt
TASK_FILE=$OUTDIR/rohan_tasks_all.tsv
ROHMU=$(cat "$OUTDIR/chosen_rohmu.txt")
CHAINS=10000

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_FILE")
IFS=$'\t' read -r prefix bam <<< "$LINE"

echo "Task ${SLURM_ARRAY_TASK_ID}: prefix=$prefix bam=$bam rohmu=$ROHMU"

if [[ -s "${prefix}.summary.txt" ]]; then
    echo "Pass 2 already complete, skipping."
    exit 0
fi

if [[ ! -s "${prefix}.hEst.gz" ]] || ! gzip -t "${prefix}.hEst.gz" 2>/dev/null; then
    echo "Pass 1 not done yet for this sample (no dependency gate — see script header) — skipping cleanly, not an error." >&2
    echo "Will be picked up on the next launch_rohan.sh resubmission once phase 1 finishes." >&2
    exit 0
fi

"$ROHAN" \
    --tstv 1.86 --size 100000 --chains "$CHAINS" \
    --hmm --rohmu "$ROHMU" \
    --auto "$SCAFFOLDS" \
    -o "$prefix" "$REF" "$bam" \
    2>> "${prefix}.log"

echo "Pass 2 complete."
