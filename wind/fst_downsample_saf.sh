#!/bin/bash
# SLURM array task, one population x replicate (submitted by fst_downsample_launch.sh).
# Builds the SAF for one n = 4 subset with the same settings as pop_SAF.sh.
# Skips a SAF that already exists.
# Environment: angsd_dev; expects OUTDIR exported by the launcher.

#SBATCH --job-name=fstds_saf
#SBATCH --partition=long
#SBATCH --cpus-per-task=9
#SBATCH --mem=32gb
#SBATCH --time=168:00:00

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

ANGSD=${BASE_DIR}/software/angsd_dev/angsd
REF=${BASE_DIR}/genome/DEINVI_genome.fasta.masked
SITES=${BASE_DIR}/selection/TajD/all_sites
CHRS=${BASE_DIR}/selection/TajD/all_chrs.txt

read -r POP REP <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${OUTDIR}/saf_tasks.txt")"
OUT="${OUTDIR}/saf/${POP}_rep${REP}"
echo "Task ${SLURM_ARRAY_TASK_ID}: ${POP}, replicate ${REP}"
cat "${OUTDIR}/bamlists/${POP}_rep${REP}.txt"

if [ -s "${OUT}.saf.idx" ]; then
  echo "SAF exists for ${POP} rep ${REP}, skipping"
  exit 0
fi

"${ANGSD}" -bam "${OUTDIR}/bamlists/${POP}_rep${REP}.txt" -ref "${REF}" -anc "${REF}" \
  -out "${OUT}" -nThreads 9 \
  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
  -minMapQ 20 -minQ 20 -doSaf 1 -GL 1 \
  -sites "${SITES}" -rf "${CHRS}"
echo "$(date) -- Done SAF: ${POP} rep ${REP}"
