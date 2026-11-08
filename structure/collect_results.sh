#!/bin/bash
# collect_results.sh
# Concatenate per-pair .stat files into pairwise_Fst.txt (submitted by
# calculate_Fst_individuals.sh).

#SBATCH --job-name=Fst_collect
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=1:00:00
#SBATCH --output=logs/collect_%j.log

BASE_DIR=/private/groups/meyerlab/esteiger

DIR="${BASE_DIR}/structure/Fst/individuals"
PAIRS_DIR="${DIR}/pairs"
OUT="${DIR}/pairwise_Fst.txt"

N_STAT=$(find "${PAIRS_DIR}" -name "*.stat" | wc -l)
echo "Collecting ${N_STAT} pair results into ${OUT}"

N_PAIRS=$(wc -l < "${DIR}/pairs.txt")
if [ "${N_STAT}" -ne "${N_PAIRS}" ]; then
    echo "WARNING: found ${N_STAT} .stat files but expected ${N_PAIRS} pairs"
fi

{
    echo -e "sample1\tsample2\tFst_unweighted\tFst_weighted"
    find "${PAIRS_DIR}" -name "*.stat" | xargs cat
} > "${OUT}"

echo "Done: ${OUT}"
