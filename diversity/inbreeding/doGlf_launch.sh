#!/bin/bash
# doGlf_launch.sh
# Index the ngsF site list and submit createGlf_ngsF.sh for each BAM-list chunk
# without a .glf.gz. Usage: bash doGlf_launch.sh
set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

ANGSD_DIR="${BASE_DIR}/variant_calling/angsd"
NGSF_DIR="${BASE_DIR}/diversity/ngsF"
SITES="${ANGSD_DIR}/subsample.sites"
CHROM="${ANGSD_DIR}/subsample.chrom"

# --- Stage 1: site list (already run; kept for reference) ---
# MAF filter is applied before the individual subset (vcftools behaviour).
# Retaining data at low-coverage samples is prioritised over LD-based filtering.
# if [ ! -f "$SITES" ]; then
#   VCF="${ANGSD_DIR}/postNGSparalog_beforeMAF_auto.vcf"
#   TO_KEEP="BASE_DIR/mapping/merged/bamlist_lowcov.txt"
#
#   vcftools --vcf "$VCF" --maf 0.05 --keep "$TO_KEEP" \
#     --recode --recode-INFO-all --out "${NGSF_DIR}/filtered_vcf"
#
#   # max-missing 0.45: genotypes in >= 45% of individuals (~530k sites retained)
#   vcftools --vcf "${NGSF_DIR}/filtered_vcf.recode.vcf" --max-missing 0.45 \
#     --recode --recode-INFO-all --out "${NGSF_DIR}/missingness"
#
#   grep -v '^#' "${NGSF_DIR}/missingness.recode.vcf" | cut -f1,2 > "$SITES"
#   cut -f1 "$SITES" > "$CHROM"
#
#   echo "Stage 1 complete: $(wc -l < "$SITES") sites retained"
# else
#   echo "Stage 1 skipped: $SITES already exists"
# fi

# --- Stage 1.5: Index sites for ANGSD ---
if [ ! -f "${SITES}.bin" ]; then
  angsd sites index "$SITES"
  echo "ANGSD sites index complete"
else
  echo "ANGSD sites index skipped: ${SITES}.bin already exists"
fi

# --- Stage 2: Launch ANGSD GLF generation on BAM-list chunks ---
cd "$NGSF_DIR"

for chunk in bamlist_chunk_*[!.]*; do
  # Skip files with extensions (e.g. .glf.gz, .mafs.gz)
  [[ "$chunk" == *.* ]] && continue
  if [ ! -f "${chunk}.glf.gz" ]; then
    sbatch ${BASE_DIR}/gvtp_popgen/diversity/inbreeding/createGlf_ngsF.sh "$chunk"
  else
    echo "Skipping $chunk: ${chunk}.glf.gz already exists"
  fi
done
