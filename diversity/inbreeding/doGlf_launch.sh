#!/bin/bash
set -euo pipefail

ANGSD_DIR="/private/groups/meyerlab/esteiger/variant_calling/angsd"
NGSF_DIR="/private/groups/meyerlab/esteiger/diversity/ngsF"
SITES="${ANGSD_DIR}/subsample.sites"
CHROM="${ANGSD_DIR}/subsample.chrom"

# --- STAGE 1: Filter VCF for sites with low missingness in low-coverage samples ---
# MAF filter applied first, then individual subset (verified with vcftools behavior).
# Prioritize retaining data at low-cov samples over LD-based filtering.
# if [ ! -f "$SITES" ]; then
#   VCF="${ANGSD_DIR}/postNGSparalog_beforeMAF_auto.vcf"
#   TO_KEEP="/private/groups/meyerlab/esteiger/mapping/merged/bamlist_lowcov.txt"
#
#   vcftools --vcf "$VCF" --maf 0.05 --keep "$TO_KEEP" \
#     --recode --recode-INFO-all --out "${NGSF_DIR}/filtered_vcf"
#
#   # max-missing 0.45 = at least 45% of individuals must have genotypes (~530k sites retained)
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

# --- STAGE 1.5: Index sites for ANGSD ---
if [ ! -f "${SITES}.bin" ]; then
  angsd sites index "$SITES"
  echo "ANGSD sites index complete"
else
  echo "ANGSD sites index skipped: ${SITES}.bin already exists"
fi

# --- STAGE 2: Launch ANGSD GLF generation on bamfile chunks ---
cd "$NGSF_DIR"

for chunk in bamlist_chunk_*[!.]*; do
  # Skip files that have extensions (e.g. .glf.gz, .mafs.gz)
  [[ "$chunk" == *.* ]] && continue
  if [ ! -f "${chunk}.glf.gz" ]; then
    sbatch ./createGlf_ngsF.sh "$chunk"
  else
    echo "Skipping $chunk: ${chunk}.glf.gz already exists"
  fi
done
