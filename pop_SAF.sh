#!/bin/bash

#SBATCH --job-name=pop_SAF
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=esteiger@ucsc.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36    # 4 parallel populations × 9 threads each
#SBATCH --mem=100gb
#SBATCH --time=168:00:00
#SBATCH --output=pop_SAF_%j.log

cd /private/groups/meyerlab/esteiger/mapping/merged

OUTDIR=/private/groups/meyerlab/esteiger/SAF/population
REF=/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked
SITES=/private/groups/meyerlab/esteiger/selection/TajD/all_sites
# all_sites includes both variable (non-MAF-filtered) and invariable sites passing quality
# filters. Including invariable sites is essential so that the SFS denominator reflects the
# true number of sites examined, avoiding upward bias in heterozygosity and Fst estimates.
# No MAF filter is applied: the ANGSD genotype-likelihood framework handles sequencing errors
# probabilistically, and MAF filtering would truncate rare variants, distorting the SFS.
CHRS=/private/groups/meyerlab/esteiger/selection/TajD/all_chrs.txt

mkdir -p "${OUTDIR}"

declare -A BAMLISTS
BAMLISTS=(
  [arguello]="bamlist_4x_arguello.txt"
  [conception]="bamlist_4x_conception.txt"
  [lionhead]="bamlist_4x_lionhead.txt"
  [cdfw]="bamlist_4x_cdfw.txt"
  # we use a subset of 16 high-coverage individuals for our SAF computation, so that Strauss is equally represented relative to the other populations
  [strauss]="bamlist_4x_strauss_top16.txt"
  [Dangermond]="bamlist_Dangermond.txt"
  # per-sampling-year cohorts within Strauss (Tranquillon Mtn), >4x-depth individuals only,
  # to compare Ne across years for the same population
  [strauss_year1]="bamlist_year1_4x_strauss.txt"
  [strauss_year2]="bamlist_year2_4x_strauss.txt"
  [strauss_year3]="bamlist_year3_4x_strauss.txt"
)
# Temporarily restricted to the new per-year cohorts so the script only computes those SAFs
# (the population-level SAFs below already exist and would be skipped anyway, but commenting
# them out avoids re-checking/re-listing them). Uncomment the line below to restore the full run.
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
