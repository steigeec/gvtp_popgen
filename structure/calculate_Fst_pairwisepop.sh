#!/bin/bash
# calculate_Fst_pairwisepop.sh
# Pairwise Fst between populations. Writes pairs.txt from POPS and submits
# pop_pairwise_array.sh (one task per pair) followed by pop_collect_results.sh.
# Safe to re-run; pairs with a completed .stat are skipped.
#
# Usage: bash calculate_Fst_pairwisepop.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

# === Paths ===
export SAFDIR="${BASE_DIR}/SAF/population"
export OUTDIR="${BASE_DIR}/structure/Fst/populations/pairwise_auto"
export RS="${BASE_DIR}/software/angsd_dev/misc/realSFS"
SCRIPTS_DIR="${BASE_DIR}/gvtp_popgen/structure"

mkdir -p "${OUTDIR}" "${OUTDIR}/logs"

# SAF files are built by gvtp_popgen/pop_SAF.sh
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
