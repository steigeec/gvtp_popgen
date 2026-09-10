#!/bin/bash
# Array task (auto-submitted by calculate_Fst_pairwisepop.sh, one pair per task):
# computes 2D-SFS, Fst index, and Fst stats for one population pair.
#
# Safe to re-run: skips entirely if this pair already has a completed .stat
# result, and skips 2D-SFS / fst-index recomputation for any pair whose .ml /
# .fst.idx already exists and is non-empty — so previously-completed work
# (e.g. the strauss vs Dangermond pair) is never overwritten, only its cheap
# read-only fst stats step is (re)run to populate the .stat file.

#SBATCH --job-name=popFst_pair
#SBATCH --partition=long
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=esteiger@ucsc.edu
#SBATCH --cpus-per-task=9
#SBATCH --mem=150gb
#SBATCH --time=168:00:00
#SBATCH --output=logs/pair_%A_%a.log
# Array size set dynamically by the launcher: --array=1-N_PAIRS

# === Paths ===
SAFDIR="/private/groups/meyerlab/esteiger/SAF/population"
OUTDIR="/private/groups/meyerlab/esteiger/structure/Fst/populations/pairwise_auto"
# Use the stable realSFS, not angsd_dev's: for the 5 largest-dimension pairs
# angsd_dev's realSFS ("Multi SFS is still under development") ran the full EM
# optimization successfully but then silently produced empty 2D-SFS output
# (confirmed not an OOM — MaxRSS was ~66GB of the 150GB budget). The stable
# binary produces correct output for the same pairs and its fst.idx has been
# verified readable by angsd_dev's realSFS, so window_Fst.sh is unaffected.
RS="/private/groups/meyerlab/esteiger/software/angsd/misc/realSFS"

read -r P1 P2 <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${OUTDIR}/pairs.txt")"
echo "Task ${SLURM_ARRAY_TASK_ID}: ${P1} vs ${P2}"

RESULT="${OUTDIR}/${P1}.${P2}.stat"
LOCK="${OUTDIR}/${P1}.${P2}.lock"
ML="${OUTDIR}/${P1}.${P2}.ml"
FST="${OUTDIR}/${P1}.${P2}"

# Already fully done — never touch it again.
if [ -s "${RESULT}" ]; then
    echo "  Already have a result for ${P1} vs ${P2}, skipping"
    exit 0
fi

# Atomic lock — avoid two tasks racing on the same pair
if ! mkdir "${LOCK}" 2>/dev/null; then
    echo "  Lock exists for ${P1} vs ${P2} (another task is working on it), exiting"
    exit 0
fi
trap 'rm -rf "'"${LOCK}"'"' EXIT

# 2D-SFS (skip only if a previous run already produced a real, non-empty one)
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

# Fst index (skip only if a previous run already produced a real, non-empty one)
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

# Fst stats (cheap, read-only — safe to (re)run even for an already-complete pair)
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
