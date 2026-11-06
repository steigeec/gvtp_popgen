#!/bin/bash
# =============================================================================
# resubmit_phase2.sh — periodically resubmit a dependency-free phase2 array
# to pick up samples whose phase 1 has finished since the last resubmission.
#
# Background: rohan_phase2.sh used to be triggered via a SLURM aftercorr
# dependency on phase1, but that proved unreliable at this array size (see
# rohan_phase2.sh / launch_rohan.sh headers) — most elements never left
# "Dependency" state even days after their phase-1 counterpart finished.
# The fix was to drop the dependency entirely: rohan_phase2.sh now
# self-validates (checks its sample's hEst.gz) and exits 0 cleanly if phase 1
# isn't done yet, so phase2 just needs to be resubmitted periodically to
# sweep up newly-finished samples. This script is that periodic sweep.
#
# SAFE to run anytime, including while phase1 is still actively running —
# unlike launch_rohan.sh (which also resubmits phase1 and can conflict with
# an in-progress phase1 run on the same samples), this script ONLY touches
# phase2. phase2 only reads hEst.gz, never writes it, so there's no
# concurrent-write risk with an active phase1 job.
#
# Uses rohan_tasks_all.tsv — ALL 423 eligible samples, unconditionally, kept
# deliberately separate from phase1's rohan_tasks.tsv (which phase1 reads BY
# LINE NUMBER at dispatch time; mutating it while phase1 still has pending
# tasks would silently corrupt future dispatches). Every task here
# self-validates and skips cleanly if its sample isn't phase1-ready or is
# already done, so re-running this against the full cohort repeatedly is
# always safe, just occasionally a no-op for a given sample.
#
# Usage: bash resubmit_phase2.sh   (intended to run on a schedule — see cron)
# =============================================================================

set -euo pipefail

OUTDIR=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan
SCRIPT_DIR=/private/groups/meyerlab/esteiger/gvtp_popgen/diversity/ROHan
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
