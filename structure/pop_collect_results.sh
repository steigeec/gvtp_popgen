#!/bin/bash
# pop_collect_results.sh
# Concatenate per-pair .stat files into pairwise_Fst_populations.txt (submitted
# by calculate_Fst_pairwisepop.sh).

#SBATCH --job-name=popFst_collect
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=1:00:00
#SBATCH --output=logs/collect_%j.log

BASE_DIR=/private/groups/meyerlab/esteiger

OUTDIR="${BASE_DIR}/structure/Fst/populations/pairwise_auto"
OUT="${OUTDIR}/pairwise_Fst_populations.txt"

N_STAT=$(find "${OUTDIR}" -maxdepth 1 -name "*.stat" | wc -l)
N_PAIRS=$(wc -l < "${OUTDIR}/pairs.txt")
echo "Collecting ${N_STAT} / ${N_PAIRS} pair results into ${OUT}"
if [ "${N_STAT}" -ne "${N_PAIRS}" ]; then
    echo "WARNING: found ${N_STAT} .stat files but expected ${N_PAIRS} pairs"
fi

{
    echo -e "pop1\tpop2\tFst_unweighted\tFst_weighted"
    find "${OUTDIR}" -maxdepth 1 -name "*.stat" | sort | xargs cat
} > "${OUT}"

echo "Done: ${OUT}"
