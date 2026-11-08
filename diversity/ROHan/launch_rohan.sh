#!/bin/bash
# launch_rohan.sh
# Build the task list of samples without a summary.txt and submit the ROHan
# arrays rohan_phase1.sh and rohan_phase2.sh, one task per sample. Re-run until
# every sample is done; the task list is rebuilt each time.
#
# Phase 1 (MCMC, --nohmm) is the expensive step. It runs with -t 8, which made
# visibly more MCMC progress than -t 1 in a benchmark on a slow sample, so
# concurrency is capped lower. The 3-day limit is well above observed runtimes
# (< 24 h) but backfills more easily than the 14-day partition cap. ROHan has
# no checkpointing, so a timed-out sample restarts on the next run.
# Phase 2 (HMM, --hmm) takes minutes.
#
# Requires, in order: prep_samples.sh (samples_gt3x.tsv), estimate_rohmu.sh
# (chosen_rohmu.txt, review before continuing).
#
# Phase 2 is submitted without a dependency because aftercorr dependencies did
# not release reliably at this array size. rohan_phase2.sh exits 0 if its
# sample's hEst.gz is not ready, and the sample is picked up on the next run.
#
# Do not re-run while a phase 1 array still has running tasks: a second task
# for the same sample would write to the same prefix. To resubmit phase 2 only,
# using the existing rohan_tasks.tsv:
#   n=$(wc -l < rohan_tasks.tsv); sbatch --array=1-${n}%40 rohan_phase2.sh
#
# Logs of completed samples are gzipped (mostly progress-bar redraws).
#
# Usage: bash launch_rohan.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

OUTDIR=${BASE_DIR}/diversity/ROH/ROHan
SCRIPT_DIR=${BASE_DIR}/gvtp_popgen/diversity/ROHan
SAMPLES_TSV=$OUTDIR/samples_gt3x.tsv
ROHMU_FILE=$OUTDIR/chosen_rohmu.txt
TASK_FILE=$OUTDIR/rohan_tasks.tsv

# Phase 1: 8 cores/task. Phase 2: single-threaded.
MAX_CONCURRENT_P1=25
MAX_CONCURRENT_P2=40

# --- Validate prerequisites --------------------------------------------------
[[ -s "$SAMPLES_TSV" ]] || {
    echo "ERROR: $SAMPLES_TSV not found. Run prep_samples.sh first." >&2; exit 1
}
[[ -s "$ROHMU_FILE" ]] || {
    echo "ERROR: $ROHMU_FILE not found. Run estimate_rohmu.sh first." >&2; exit 1
}

ROHMU=$(cat "$ROHMU_FILE")
echo "rohmu: $ROHMU  (from chosen_rohmu.txt)"

# --- Build remaining-sample task list ----------------------------------------
# samples_gt3x.tsv columns: sample, cov, orig_bam, final_bam (written by prep_samples.sh)
> "$TASK_FILE"
n_total=0
n_done=0
n_compressed=0
while IFS=$'\t' read -r sample cov orig_bam final_bam; do
    n_total=$((n_total + 1))
    prefix="$OUTDIR/$sample"
    if [[ -s "${prefix}.summary.txt" ]]; then
        n_done=$((n_done + 1))
        if [[ -f "${prefix}.log" ]]; then
            gzip "${prefix}.log"
            n_compressed=$((n_compressed + 1))
        fi
        continue
    fi
    printf '%s\t%s\n' "$prefix" "$final_bam" >> "$TASK_FILE"
done < "$SAMPLES_TSV"

n_remaining=$(wc -l < "$TASK_FILE")
echo "Samples total: $n_total | already complete: $n_done | remaining: $n_remaining"
(( n_compressed > 0 )) && echo "Compressed ${n_compressed} newly-completed sample log(s) (.log -> .log.gz)"

if (( n_remaining == 0 )); then
    echo "All samples already have a summary.txt — nothing to submit."
    echo "Next: bash collect_rohan.sh, then bash plot_rohan.sh"
    exit 0
fi

mkdir -p "$OUTDIR/logs"

# --- Submit two-phase job arrays ----------------------------------------------
echo ""
echo "Submitting phase 1 (MCMC) array: $n_remaining tasks, ${MAX_CONCURRENT_P1} concurrent max (8 cores each)"
P1=$(sbatch --parsable --array=1-${n_remaining}%${MAX_CONCURRENT_P1} "$SCRIPT_DIR/rohan_phase1.sh")
echo "  Phase 1 job ID: $P1"

echo "Submitting phase 2 (HMM) array: no dependency (self-validates hEst.gz — see header note), ${MAX_CONCURRENT_P2} concurrent max"
P2=$(sbatch --parsable --array=1-${n_remaining}%${MAX_CONCURRENT_P2} "$SCRIPT_DIR/rohan_phase2.sh")
echo "  Phase 2 job ID: $P2"

echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Re-run this script at any time to pick up samples that timed out or"
echo "haven't started yet — rohan_tasks.tsv is rebuilt fresh each time."
echo "Once all 423 samples have a summary.txt: bash collect_rohan.sh, then bash plot_rohan.sh"
