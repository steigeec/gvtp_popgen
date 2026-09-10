#!/bin/bash
# =============================================================================
# rohan_phase2.sh — SLURM array: ROHan main-run pass 2 (HMM ROH calling)
#
# One array task per sample, same indexing as rohan_phase1.sh (task list is
# rohan_tasks.tsv).  Submitted with --dependency=aftercorr:<phase1 jobid> so
# each task only waits on its OWN phase-1 counterpart, not the whole phase-1
# array — see launch_rohan.sh.  Reads the empirically calibrated rohmu value
# from chosen_rohmu.txt (written by estimate_rohmu.sh).
#
# Fast relative to phase 1: HMM calling with a single fixed rohmu value took
# well under 10 minutes per sample during rohmu calibration (chains=5000);
# the --time budget below has generous headroom over that.
# =============================================================================

#SBATCH --job-name=rohan_p2
#SBATCH --partition=medium
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=esteiger@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=06:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/logs/rohan_phase2_%A_%a.log

set -euo pipefail

ROHAN=/private/groups/meyerlab/esteiger/software/rohan/bin/rohan
REF=/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked
OUTDIR=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan
SCAFFOLDS=$OUTDIR/scaffolds_1mb_names.txt
TASK_FILE=$OUTDIR/rohan_tasks.tsv
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
    echo "ERROR: ${prefix}.hEst.gz missing or invalid — pass 1 did not complete for this sample." >&2
    echo "This task should not have run (aftercorr dependency); skipping rather than feeding the HMM garbage." >&2
    exit 1
fi

"$ROHAN" \
    --tstv 1.86 --size 100000 --chains "$CHAINS" \
    --hmm --rohmu "$ROHMU" \
    --auto "$SCAFFOLDS" \
    -o "$prefix" "$REF" "$bam" \
    2>> "${prefix}.log"

echo "Pass 2 complete."
