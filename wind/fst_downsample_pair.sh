#!/bin/bash
# SLURM array task, one pair x replicate (submitted by fst_downsample_launch.sh).
# Folded 2D-SFS and Fst for one population pair, as in structure/pop_pairwise_array.sh.
# Dangermond uses its full-sample SAF (already n = 4).
# Writes pairs/<pop1>.<pop2>.rep<r>.stat (pop1, pop2, replicate, unweighted Fst, weighted Fst).
# Environment: stable angsd realSFS; expects OUTDIR and FULL_SAFDIR exported by the launcher.

#SBATCH --job-name=fstds_pair
#SBATCH --partition=long
#SBATCH --cpus-per-task=9
#SBATCH --mem=64gb
#SBATCH --time=72:00:00

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

# Stable realSFS; the angsd_dev build gave empty 2D-SFS output for some pairs.
RS=${BASE_DIR}/software/angsd/misc/realSFS

read -r P1 P2 REP <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${OUTDIR}/pair_tasks.txt")"
echo "Task ${SLURM_ARRAY_TASK_ID}: ${P1} vs ${P2}, replicate ${REP}"
saf() { if [ "$1" = Dangermond ]; then echo "${FULL_SAFDIR}/Dangermond.saf.idx"; else echo "${OUTDIR}/saf/$1_rep${REP}.saf.idx"; fi; }

BASE="${OUTDIR}/pairs/${P1}.${P2}.rep${REP}"
if [ -s "${BASE}.stat" ]; then
  echo "Already done, skipping"
  exit 0
fi

"${RS}" "$(saf "${P1}")" "$(saf "${P2}")" -P 9 -fold 1 > "${BASE}.ml"
[ -s "${BASE}.ml" ] || { echo "ERROR: 2D-SFS failed"; exit 1; }
"${RS}" fst index "$(saf "${P1}")" "$(saf "${P2}")" -sfs "${BASE}.ml" -fstout "${BASE}"
STATS=$("${RS}" fst stats "${BASE}.fst.idx" | grep -v '^\s*->' | grep -v '^$' | tail -1)
UNWT=$(echo "${STATS}" | awk '{print $1}')
WGT=$(echo "${STATS}" | awk '{print $2}')
if [[ "${UNWT}" =~ ^-?[0-9.eE+-]+$ ]] && [[ "${WGT}" =~ ^-?[0-9.eE+-]+$ ]]; then
  echo -e "${P1}\t${P2}\t${REP}\t${UNWT}\t${WGT}" > "${BASE}.stat"
  rm -f "${BASE}.fst.gz" "${BASE}.fst.idx"   # per-site output is large
  echo "Done: Fst unweighted ${UNWT}, weighted ${WGT}"
else
  echo "ERROR: Fst stats failed"
  exit 1
fi
