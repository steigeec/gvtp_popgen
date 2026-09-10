#!/bin/bash
# Launcher: generates pairs.txt from POPS, then submits one array task per pair
# (150gb/9cpus each, scheduled independently instead of one big single-node job)
# followed by a collect job as a SLURM dependency chain.
#
# Safe to re-run: pop_pairwise_array.sh skips any pair that already has a
# completed .stat result, so previously-finished pairs (e.g. strauss vs
# Dangermond) are never recomputed or overwritten.
#
# Usage: bash calculate_Fst_pairwisepop.sh

set -euo pipefail

# === Paths ===
export SAFDIR="/private/groups/meyerlab/esteiger/SAF/population"
export OUTDIR="/private/groups/meyerlab/esteiger/structure/Fst/populations/pairwise_auto"
export RS="/private/groups/meyerlab/esteiger/software/angsd_dev/misc/realSFS"
SCRIPTS_DIR="/private/groups/meyerlab/esteiger/gvtp_popgen/structure"

mkdir -p "${OUTDIR}" "${OUTDIR}/logs"

# SAF files are built by /private/groups/meyerlab/esteiger/gvtp_popgen/pop_SAF.sh
POPS=(conception cdfw arguello strauss lionhead Dangermond)

# === Generate all unique pairs (deterministic) ===
: > "${OUTDIR}/pairs.txt"
for (( i=0; i<${#POPS[@]}-1; i++ )); do
  for (( j=i+1; j<${#POPS[@]}; j++ )); do
    echo "${POPS[$i]} ${POPS[$j]}" >> "${OUTDIR}/pairs.txt"
  done
done
N_PAIRS=$(wc -l < "${OUTDIR}/pairs.txt")
echo "$(date) -- ${N_PAIRS} pairwise comparisons"

# === Submit one array task per pair ===
ARR_JID=$(sbatch \
    --job-name=popFst_pair \
    --partition=long \
    --mail-type=FAIL \
    --mail-user=esteiger@ucsc.edu \
    --array=1-${N_PAIRS} \
    --cpus-per-task=9 \
    --mem=150gb \
    --time=168:00:00 \
    --output="${OUTDIR}/logs/pair_%A_%a.log" \
    --export=ALL \
    "${SCRIPTS_DIR}/pop_pairwise_array.sh" | awk '{print $NF}')
echo "Submitted pairwise array job: ${ARR_JID} (${N_PAIRS} tasks, 150gb/9cpus each)"

# === Submit collect job (waits for all array tasks) ===
COL_JID=$(sbatch \
    --job-name=popFst_collect \
    --partition=short \
    --mail-type=END,FAIL \
    --mail-user=esteiger@ucsc.edu \
    --ntasks=1 \
    --mem=4gb \
    --time=1:00:00 \
    --output="${OUTDIR}/logs/collect_%j.log" \
    --dependency=afterok:${ARR_JID} \
    --export=ALL \
    "${SCRIPTS_DIR}/pop_collect_results.sh" | awk '{print $NF}')
echo "Submitted collect job: ${COL_JID} (depends on ${ARR_JID})"
echo ""
echo "Pipeline: ${ARR_JID} → ${COL_JID}"
echo "Final table will be written to: ${OUTDIR}/pairwise_Fst_populations.txt"
