#!/bin/bash

#SBATCH --job-name=pop_SAF
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36    # 4 parallel populations × 9 threads each
#SBATCH --mem=100gb
#SBATCH --time=168:00:00
#SBATCH --output=pop_SAF_%j.log

# Per-population SAFs (ANGSD -doSaf) over all_sites, for SFS, heterozygosity,
# Fst and Watterson's theta.

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/mapping/merged

OUTDIR=${BASE_DIR}/SAF/population
REF=${BASE_DIR}/genome/DEINVI_genome.fasta.masked
SITES=${BASE_DIR}/selection/TajD/all_sites
# all_sites includes both variable (no MAF filter) and invariable sites passing quality
# filters. Invariable sites are needed so the SFS denominator reflects the true number of
# sites examined, avoiding upward bias in heterozygosity and Fst estimates.
# No MAF filter: ANGSD's genotype-likelihood framework handles sequencing errors
# probabilistically, and a MAF filter would truncate rare variants and distort the SFS.
CHRS=${BASE_DIR}/selection/TajD/all_chrs.txt

mkdir -p "${OUTDIR}"

declare -A BAMLISTS
BAMLISTS=(
  [arguello]="bamlist_4x_arguello.txt"
  [conception]="bamlist_4x_conception.txt"
  [lionhead]="bamlist_4x_lionhead.txt"
  [cdfw]="bamlist_4x_cdfw.txt"
  # 16 high-coverage individuals, so Strauss is represented equally with the other populations
  [strauss]="bamlist_4x_strauss_top16.txt"
  [Dangermond]="bamlist_Dangermond.txt"
  # per-sampling-year cohorts within Strauss, >4x individuals only
  [strauss_year1]="bamlist_year1_4x_strauss.txt"
  [strauss_year2]="bamlist_year2_4x_strauss.txt"
  [strauss_year3]="bamlist_year3_4x_strauss.txt"
)
# for the population-level SAFs, use this line instead of the per-year cohorts
# POPS=(arguello conception lionhead cdfw strauss Dangermond)
POPS=(strauss_year1 strauss_year2 strauss_year3)

declare -p BAMLISTS > /tmp/bamlists_export.sh
export OUTDIR REF SITES CHRS

printf '%s\n' "${POPS[@]}" | parallel -j 4 '
  POP={}
  source /tmp/bamlists_export.sh
  BAMLIST="${BAMLISTS[$POP]}"
  OUT="${OUTDIR}/${POP}"

  if [ -f "${OUT}.saf.idx" ]; then
    echo "SAF exists for ${POP}, skipping"
    exit 0
  fi

  echo "$(date) -- Computing SAF: ${POP}"
  angsd -bam "${BAMLIST}" -ref "${REF}" -anc "${REF}" \
    -out "${OUT}" -nThreads 9 \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
    -minMapQ 20 -minQ 20 -doSaf 1 -GL 1 \
    -sites "${SITES}" -rf "${CHRS}" \
    > "${OUT}.angsd.log" 2>&1
  echo "$(date) -- Done SAF: ${POP}"
'
