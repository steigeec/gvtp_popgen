#!/bin/bash
# pop_pairwise_array.sh
# realSFS 2D-SFS, Fst index and Fst stats for one population pair (submitted by
# calculate_Fst_pairwisepop.sh). Existing non-empty .stat, .ml and .fst.idx files
# are reused, not overwritten.

#SBATCH --job-name=popFst_pair
#SBATCH --partition=long
#SBATCH --cpus-per-task=9
#SBATCH --mem=150gb
#SBATCH --time=168:00:00
#SBATCH --output=logs/pair_%A_%a.log
# --array is set by the launcher.

# === Paths ===
BASE_DIR=/private/groups/meyerlab/esteiger

SAFDIR="${BASE_DIR}/SAF/population"
OUTDIR="${BASE_DIR}/structure/Fst/populations/pairwise_auto"
# Stable realSFS rather than angsd_dev's: for the 5 largest pairs angsd_dev's
# realSFS ran the EM but wrote an empty 2D-SFS (not an OOM). The stable binary's
# fst.idx is readable by angsd_dev's realSFS, so window_Fst.sh is unaffected.
RS="${BASE_DIR}/software/angsd/misc/realSFS"

read -r P1 P2 <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${OUTDIR}/pairs.txt")"
echo "Task ${SLURM_ARRAY_TASK_ID}: ${P1} vs ${P2}"

RESULT="${OUTDIR}/${P1}.${P2}.stat"
LOCK="${OUTDIR}/${P1}.${P2}.lock"
ML="${OUTDIR}/${P1}.${P2}.ml"
FST="${OUTDIR}/${P1}.${P2}"

if [ -s "${RESULT}" ]; then
    echo "  Already have a result for ${P1} vs ${P2}, skipping"
    exit 0
fi

# Atomic lock to avoid two tasks racing on the same pair
if ! mkdir "${LOCK}" 2>/dev/null; then
    echo "  Lock exists for ${P1} vs ${P2} (another task is working on it), exiting"
    exit 0
fi
trap 'rm -rf "'"${LOCK}"'"' EXIT

# 2D-SFS (skip if a non-empty one exists)
if [ ! -s "${ML}" ]; then
    rm -f "${ML}"
    echo "  2D-SFS: ${P1} vs ${P2}"
    "${RS}" "${SAFDIR}/${P1}.saf.idx" "${SAFDIR}/${P2}.saf.idx" -P 9 -fold 1 > "${ML}"
    if [ ! -s "${ML}" ]; then
        echo "  ERROR: 2D-SFS failed for ${P1} vs ${P2}"
        rm -f "${ML}"
        exit 1
    fi
else
    echo "  2D-SFS already exists for ${P1} vs ${P2}, skipping"
fi

# Fst index (skip if a non-empty one exists)
if [ ! -s "${FST}.fst.idx" ]; then
    echo "  Fst index: ${P1} vs ${P2}"
    "${RS}" fst index "${SAFDIR}/${P1}.saf.idx" "${SAFDIR}/${P2}.saf.idx" \
        -sfs "${ML}" -fstout "${FST}"
else
    echo "  Fst index already exists for ${P1} vs ${P2}, skipping"
fi

if [ ! -s "${FST}.fst.idx" ]; then
    echo "  ERROR: Fst index missing/empty for ${P1} vs ${P2}"
    exit 1
fi

# Fst stats (read-only, safe to rerun)
STATS=$("${RS}" fst stats "${FST}.fst.idx" | grep -v '^\s*->' | grep -v '^$' | tail -1)
UNWT=$(echo "${STATS}" | awk '{print $1}')
WGT=$(echo "${STATS}"  | awk '{print $2}')

if [[ "${UNWT}" =~ ^-?[0-9.eE+-]+$ ]] && [[ "${WGT}" =~ ^-?[0-9.eE+-]+$ ]]; then
    echo -e "${P1}\t${P2}\t${UNWT}\t${WGT}" > "${RESULT}.part"
    mv "${RESULT}.part" "${RESULT}"
    echo "  Done: ${P1} vs ${P2}"
else
    echo "  ERROR: Fst stats failed for ${P1} vs ${P2}, will retry next run"
    exit 1
fi
