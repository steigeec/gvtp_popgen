#!/bin/bash
# pairwise_array.sh
# realSFS 2D-SFS and unweighted/weighted Fst for each pair of individuals in this
# task's batch of pairs.txt (submitted by calculate_Fst_individuals.sh).

#SBATCH --job-name=Fst_pairwise
#SBATCH --partition=long
#SBATCH --cpus-per-task=2
# The He SAFs cover all_sites (~455M sites incl. invariant), for which realSFS
# needs ~24 GB per pair even though the per-pair SFS is only 3x3.
#SBATCH --mem=40gb
#SBATCH --time=168:00:00
#SBATCH --output=logs/pairwise_%A_%a.log
# --array is set by the launcher.

# === Paths ===
BASE_DIR=/private/groups/meyerlab/esteiger

RS="${BASE_DIR}/software/angsd/misc/realSFS"
DIR="${BASE_DIR}/structure/Fst/individuals"
BD="${BASE_DIR}/mapping/merged"
BAMLIST="${BD}/bamlist_4x.txt"
# Per-individual SAFs from the He pipeline (launch_He.sh)
SAF_DIR="${BASE_DIR}/diversity/He/individual"
PAIRS_DIR="${DIR}/pairs"

# BATCH_SIZE must match the value in calculate_Fst_individuals.sh
BATCH_SIZE=50

# === Map each 4-char sample code to its He SAF basename (bam name minus .bam) ===
declare -A SAF_BASENAME
while IFS= read -r bam; do
    [ -z "${bam}" ] && continue
    base=$(basename "${bam}")
    base="${base%.bam}"
    SAF_BASENAME["${base:0:4}"]="${base}"
done < "${BAMLIST}"

# === Compute this task's line range in pairs.txt ===
START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))
echo "Task ${SLURM_ARRAY_TASK_ID}: pairs ${START}–${END}"

# === Process each pair in this task's batch ===
while IFS=' ' read -r S1 S2; do

    RESULT="${PAIRS_DIR}/${S1}.${S2}.stat"
    LOCK="${PAIRS_DIR}/${S1}.${S2}.lock"
    TMP="${PAIRS_DIR}/${S1}.${S2}.tmp"

    if [ -f "${RESULT}" ]; then
        continue
    fi

    SAF1="${SAF_DIR}/${SAF_BASENAME[$S1]}.saf.idx"
    SAF2="${SAF_DIR}/${SAF_BASENAME[$S2]}.saf.idx"
    # A non-empty .saf.idx may still be mid-write; ind_het.sh writes <sample>.ml
    # last, so require it as the completion marker.
    ML1="${SAF_DIR}/${SAF_BASENAME[$S1]}.ml"
    ML2="${SAF_DIR}/${SAF_BASENAME[$S2]}.ml"
    if [ -z "${SAF_BASENAME[$S1]:-}" ] || [ -z "${SAF_BASENAME[$S2]:-}" ] || \
       [ ! -s "${SAF1}" ] || [ ! -s "${SAF2}" ] || \
       [ ! -s "${ML1}" ] || [ ! -s "${ML2}" ]; then
        echo "  Skipping ${S1} vs ${S2}: SAF not ready yet (${SAF1} / ${SAF2})"
        continue
    fi

    # Atomic lock; skip if another task is working on this pair
    if ! mkdir "${LOCK}" 2>/dev/null; then
        echo "  Skipping ${S1} vs ${S2} (lock exists)"
        continue
    fi

    trap "rm -rf '${LOCK}' '${TMP}.ml' '${TMP}.fst.idx' '${TMP}.fst.gz'" EXIT

    echo "  ${S1} vs ${S2}"

    # 2D SFS
    ${RS} "${SAF1}" "${SAF2}" \
        -P 2 -fold 1 > "${TMP}.ml" 2>/dev/null

    # Fst index
    ${RS} fst index \
        "${SAF1}" "${SAF2}" \
        -sfs "${TMP}.ml" -fstout "${TMP}" 2>/dev/null

    # fst stats prints "unweighted\tweighted" after the "->" lines
    STATS=$(${RS} fst stats "${TMP}.fst.idx" 2>/dev/null | grep -v '^\s*->' | grep -v '^$' | tail -1)
    UNWT=$(echo "${STATS}" | awk '{print $1}')
    WGT=$(echo "${STATS}"  | awk '{print $2}')

    # Only write a result if realSFS returned numbers, so failed pairs are retried
    if [[ "${UNWT}" =~ ^-?[0-9.eE+-]+$ ]] && [[ "${WGT}" =~ ^-?[0-9.eE+-]+$ ]]; then
        # mv is atomic within the same filesystem
        echo -e "${S1}\t${S2}\t${UNWT}\t${WGT}" > "${RESULT}.part"
        mv "${RESULT}.part" "${RESULT}"
    else
        echo "  ERROR: Fst stats failed for ${S1} vs ${S2}, will retry next run"
    fi

    rm -f "${TMP}.ml" "${TMP}.fst.idx" "${TMP}.fst.gz"
    rm -rf "${LOCK}"

    # Reset trap after manual cleanup
    trap - EXIT

done < <(sed -n "${START},${END}p" "${DIR}/pairs.txt")

echo "Task ${SLURM_ARRAY_TASK_ID} complete"
