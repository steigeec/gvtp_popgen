#!/bin/bash
# =============================================================================
# launch_rohan.sh — Build remaining-sample task list, submit the two-phase
#                    ROHan SLURM arrays (rohan_phase1.sh, rohan_phase2.sh)
#
# Prerequisites (must be run first, in order):
#   1. prep_samples.sh    → writes samples_gt3x.tsv / bamlist_rohan.txt
#   2. estimate_rohmu.sh  → writes chosen_rohmu.txt  (REVIEW before continuing)
#
# One SLURM array task per remaining sample (not per chunk) so each sample's
# progress and time limit are independent, and SLURM can backfill freely.
# Rebuilds rohan_tasks.tsv from scratch every run, skipping any sample that
# already has a summary.txt — safe to re-run this after a partial/timed-out
# pass to pick up only what's left; array indices are reassigned each time.
#
# Phase 1 (rohan_phase1.sh): MCMC het estimation (--nohmm), the expensive step
#   — up to 3 days per sample (submitted on the `long` partition, whose real
#   cap is 14 days, but observed runtimes top out well under 24h even for the
#   highest-coverage samples, so a shorter ask backfills into the scheduler
#   far more easily than a 14-day request would), MAX_CONCURRENT_P1 running at
#   once.  Runs with -t 8 (8 threads/sample): a benchmark on a representative
#   slow sample showed threads=8 reaching visibly more MCMC progress than
#   threads=1 in the same wall-clock window, so concurrency here is capped
#   lower (8 cores/task) to keep the cluster footprint reasonable.  A sample
#   that genuinely needs longer than 3 days will TIMEOUT and restart from
#   scratch on the next re-run of this script (ROHan has no checkpointing).
# Phase 2 (rohan_phase2.sh): HMM ROH calling (--hmm), fast (~minutes). Each
#   phase 2 array task depends only on its OWN phase 1 counterpart via SLURM's
#   `aftercorr` dependency — NOT on the whole phase 1 array finishing — so
#   samples that finish phase 1 early get their phase 2 immediately rather
#   than waiting on the slowest straggler.  A phase 2 task whose phase 1
#   counterpart times out or fails simply never runs (DependencyNeverSatisfied)
#   and gets cancelled; that sample is picked up again on the next re-run of
#   this script once phase 1 succeeds for it.
#
# Also gzips each completed sample's ${prefix}.log the first time it's seen
# done — nothing reads that file's content, and it's almost entirely repeated
# ASCII progress-bar redraws (no terminal to overwrite them in place), so it
# compresses ~99%.  Only touches samples with a summary.txt already, so it
# never races a still-running phase 1/2 task still appending to its log.
#
# Usage: bash launch_rohan.sh   (safe to re-run repeatedly until all done)
# =============================================================================

set -euo pipefail

OUTDIR=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan
SCRIPT_DIR=/private/groups/meyerlab/esteiger/gvtp_popgen/diversity/ROHan
SAMPLES_TSV=$OUTDIR/samples_gt3x.tsv
ROHMU_FILE=$OUTDIR/chosen_rohmu.txt
TASK_FILE=$OUTDIR/rohan_tasks.tsv

# Phase 1 runs -t 8 (8 cores/task): 25 concurrent = 200 cores footprint.
# Phase 2 is single-threaded and fast: 40 concurrent = 40 cores footprint.
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

echo "Submitting phase 2 (HMM) array: aftercorr on phase 1 (per-sample, not array-wide), ${MAX_CONCURRENT_P2} concurrent max"
P2=$(sbatch --parsable --dependency=aftercorr:"$P1" --array=1-${n_remaining}%${MAX_CONCURRENT_P2} "$SCRIPT_DIR/rohan_phase2.sh")
echo "  Phase 2 job ID: $P2"

echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Re-run this script at any time to pick up samples that timed out or"
echo "haven't started yet — rohan_tasks.tsv is rebuilt fresh each time."
echo "Once all 423 samples have a summary.txt: bash collect_rohan.sh, then bash plot_rohan.sh"
