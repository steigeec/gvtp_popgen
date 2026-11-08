#!/bin/bash
# calculate_Fst_individuals.sh
# Pairwise Fst between individuals (step 1, entry point). Writes samples.txt and
# pairs.txt, checks that every sample's SAF is ready, then submits
# pairwise_array.sh followed by collect_results.sh (afterok). Safe to re-run;
# array tasks skip finished or locked pairs.
#
# Per-individual SAFs are reused from diversity/He/individual (same BAMs), made
# by diversity/heterozygosity/individual/launch_He.sh, which must finish first.
#
# Usage: bash calculate_Fst_individuals.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

# === Paths (shared with array scripts) ===
export DIR="${BASE_DIR}/structure/Fst/individuals"
export BD="${BASE_DIR}/mapping/merged"
export BAMLIST="${BD}/bamlist_4x.txt"
export SAF_DIR="${BASE_DIR}/diversity/He/individual"
export PAIRS_DIR="${DIR}/pairs"
export SCRIPTS_DIR="${BASE_DIR}/gvtp_popgen/structure"

BATCH_SIZE=50   # pairs per pairwise array task (keep in sync with pairwise_array.sh)

mkdir -p "${PAIRS_DIR}" "${DIR}/logs"

# === Generate samples.txt ===
echo "Generating samples.txt..."
awk 'NF { print substr($1, 1, 4) }' "${BAMLIST}" > "${DIR}/samples.txt"
N=$(wc -l < "${DIR}/samples.txt")
echo "  ${N} samples"

# === Generate pairs.txt ===
echo "Generating pairs.txt..."
mapfile -t SAMPLES < "${DIR}/samples.txt"
{
    for (( i=0; i<${#SAMPLES[@]}-1; i++ )); do
        for (( j=i+1; j<${#SAMPLES[@]}; j++ )); do
            echo "${SAMPLES[$i]} ${SAMPLES[$j]}"
        done
    done
} > "${DIR}/pairs.txt"
N_PAIRS=$(wc -l < "${DIR}/pairs.txt")
N_TASKS=$(( (N_PAIRS + BATCH_SIZE - 1) / BATCH_SIZE ))
echo "  ${N_PAIRS} pairs → ${N_TASKS} pairwise array tasks (batch size ${BATCH_SIZE})"

# === Check that every sample's SAF is ready in the He pipeline's output ===
echo "Checking SAF availability in ${SAF_DIR}..."
MISSING=0
while IFS= read -r bam; do
    [ -z "${bam}" ] && continue
    base=$(basename "${bam}")
    base="${base%.bam}"
    # A non-empty .saf.idx may still be mid-write; ind_het.sh writes <sample>.ml
    # last, so require it as the completion marker.
    if [ ! -s "${SAF_DIR}/${base}.saf.idx" ] || [ ! -s "${SAF_DIR}/${base}.ml" ]; then
        MISSING=$((MISSING + 1))
    fi
done < "${BAMLIST}"

if [ "${MISSING}" -gt 0 ]; then
    echo "  ${MISSING} / ${N} samples still missing SAF — launch_He.sh has not finished for all of them yet."
    echo "  Not submitting the pairwise job. Re-run this script once launch_He.sh completes."
    exit 1
fi
echo "  All ${N} samples have SAF ready."

# === Submit pairwise array job ===
PW_JID=$(sbatch \
    --job-name=Fst_pairwise \
    --partition=long \
    --array=1-${N_TASKS}%150 \
    --cpus-per-task=2 \
    --mem=40gb \
    --time=168:00:00 \
    --output="${DIR}/logs/pairwise_%A_%a.log" \
    --export=ALL \
    "${SCRIPTS_DIR}/pairwise_array.sh" | awk '{print $NF}')
echo "Submitted pairwise job array: ${PW_JID} (${N_TASKS} tasks)"

# === Submit collect job (waits for all pairwise tasks) ===
COL_JID=$(sbatch \
    --job-name=Fst_collect \
    --partition=short \
    --ntasks=1 \
    --mem=4gb \
    --time=1:00:00 \
    --output="${DIR}/logs/collect_%j.log" \
    --dependency=afterok:${PW_JID} \
    --export=ALL \
    "${SCRIPTS_DIR}/collect_results.sh" | awk '{print $NF}')
echo "Submitted collect job: ${COL_JID} (depends on ${PW_JID})"
echo ""
echo "Pipeline: ${PW_JID} → ${COL_JID}"
echo "Final table will be written to: ${DIR}/pairwise_Fst.txt"
