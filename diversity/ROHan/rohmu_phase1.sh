#!/bin/bash
# =============================================================================
# rohmu_phase1.sh — SLURM array: ROHan pass 1 (MCMC het estimation)
#
# Each array task runs a single pass-1 job (one sample at one depth).
# Task list is read from rohmu_calibration/phase1_tasks.tsv, written by
# estimate_rohmu.sh.  Submit via estimate_rohmu.sh or manually:
#   sbatch --array=1-N rohmu_phase1.sh
# =============================================================================

#SBATCH --job-name=rohmu_p1
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=168:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/rohmu_calibration/logs/phase1_%A_%a.log

set -euo pipefail

ROHAN=/private/groups/meyerlab/esteiger/software/rohan/bin/rohan
REF=/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked
SCAFFOLDS=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/scaffolds_1mb_names.txt
TASK_FILE=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/rohmu_calibration/phase1_tasks.tsv

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
