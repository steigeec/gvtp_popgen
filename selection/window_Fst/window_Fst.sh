#!/bin/bash

#SBATCH --job-name=window_Fst
#SBATCH --partition=long
#SBATCH --cpus-per-task=27
#SBATCH --mem=100gb
#SBATCH --time=168:00:00
#SBATCH --output=window_Fst_%j.log

# === Paths ===
# .ml and .fst.idx files are built by /private/groups/meyerlab/esteiger/gvtp_popgen/structure/calculate_Fst_pairwisepop.sh
SCRIPT_DIR="/private/groups/meyerlab/esteiger/gvtp_popgen/selection/window_Fst"
FSTDIR="/private/groups/meyerlab/esteiger/structure/Fst/populations/pairwise_auto"
PAIR_TABLE="${FSTDIR}/pairwise_Fst_populations.txt"
OUTDIR="/private/groups/meyerlab/esteiger/selection/window_Fst"
RS="/private/groups/meyerlab/esteiger/software/angsd_dev/misc/realSFS"

mkdir -p "${OUTDIR}"

export FSTDIR OUTDIR RS

# === Read all pairs from the pairwise Fst table (skip header) ===
PAIRS=()
while IFS=$'\t' read -r P1 P2 _ _; do
  PAIRS+=("${P1} ${P2}")
done < <(tail -n +2 "${PAIR_TABLE}")
echo "$(date) -- ${#PAIRS[@]} pairwise comparisons"

# === Sliding window Fst, per pair (3 at a time, 9 threads each = 27 cores) ===
echo "$(date) -- Sliding window Fst (win 50000, step 10000)"

printf '%s\n' "${PAIRS[@]}" | parallel -j 3 '
  P1=$(echo {} | cut -d" " -f1)
  P2=$(echo {} | cut -d" " -f2)
  FST="${FSTDIR}/${P1}.${P2}.fst.idx"
  OUT="${OUTDIR}/${P1}.${P2}_window"

  if [ ! -s "${FST}" ]; then
    echo "  ERROR: missing fst index for ${P1} vs ${P2}: ${FST}" >&2
    exit 1
  fi

  echo "  window Fst: ${P1} vs ${P2}" >&2
  ${RS} fst stats2 "${FST}" -win 50000 -step 10000 -P 9 -fold 1 > "${OUT}"
'

# === Plot each pairwise comparison as a genome-wide Fst scatter ===
echo "$(date) -- Plotting window Fst"
Rscript "${SCRIPT_DIR}/plot_window_Fst.R" "${OUTDIR}" "${PAIR_TABLE}"

echo "$(date) -- Done! Results in ${OUTDIR}"
