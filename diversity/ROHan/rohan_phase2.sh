#!/bin/bash
# =============================================================================
# rohan_phase2.sh — SLURM array: ROHan main-run pass 2 (HMM ROH calling)
#
# One array task per sample. Task list is rohan_tasks_all.tsv — ALL 423
# eligible samples, unconditionally, built once from samples_gt3x.tsv.
# This is DELIBERATELY a separate, stable file from rohan_phase1.sh's
# rohan_tasks.tsv (which only lists whatever subset was still incomplete
# when launch_rohan.sh last rebuilt it). Phase 1 reads its task file BY LINE
# NUMBER at dispatch time, so mutating that file while a phase-1 array still
# has pending (not-yet-started) tasks would silently corrupt future
# dispatches (wrong sample per index) — confirmed this matters in practice
# 2026-09-10 when reprocessing everyone at a new rohmu while phase1 was
# still mid-run. Phase 2 has no such constraint (each task is independent,
# self-validated below), so it gets its own full-cohort, rarely-changing
# task file instead of sharing phase1's.
#
# Submitted with NO dependency — an earlier version used
# --dependency=aftercorr:<phase1 jobid> to wait on each task's own phase-1
# counterpart, but aftercorr proved unreliable at this array size in this
# SLURM environment (most elements never transitioned out of "Dependency"
# even after their phase-1 counterpart completed days earlier — confirmed
# via `scontrol show job`). Instead this script self-validates (below) and
# exits cleanly if its sample's phase 1 isn't done yet or is already done;
# resubmit_phase2.sh re-submits a fresh array (no dependency) periodically,
# so stragglers and newly-finished phase1 samples just get picked up on the
# next sweep. Reads the calibrated rohmu value from chosen_rohmu.txt.
#
# Resource sizing (checked empirically, 2026-09-09, against 466 real HMM
# runs from the phase2 backlog-catchup batch, rohmu=1e-5): max elapsed
# 17m18s, max MaxRSS ~10MB. --mem still gives huge headroom over that.
# --time bumped after rohmu=2e-3 reprocessing (2026-09-10) hit 7 timeouts at
# the previous 45min budget — most samples still finish in ~14min regardless
# of rohmu, but a sample-specific minority apparently needs meaningfully
# longer at the higher rohmu (mechanism not yet investigated). partition
# bumped to medium (12h cap) to allow the longer --time; short (1h cap) no
# longer fits.
# =============================================================================

#SBATCH --job-name=rohan_p2
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=512mb
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/logs/rohan_phase2_%A_%a.log

set -euo pipefail

ROHAN=/private/groups/meyerlab/esteiger/software/rohan/bin/rohan
REF=/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked
OUTDIR=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan
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
