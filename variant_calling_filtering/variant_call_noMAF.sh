#!/bin/bash
# ANGSD SNP calling (GLs, mafs, beagle, BCF, genotypes) for one genome chunk,
# with no MAF filter.
# Usage: variant_call_noMAF.sh <chunk_regions_file> <bamlist> <ref.fasta>  (run by launch_variant_call_noMAF.sh)

BASE_DIR=/private/groups/meyerlab/esteiger

chunk=$1
BAMS=$2
REF=$3
indir="${BASE_DIR}/variant_calling/angsd/chromo_30_chunks/"
outdir="${BASE_DIR}/variant_calling/angsd"
export TMPDIR="/tmp"
prefix="./"
string=${chunk#"$prefix"}
ANGSD="${BASE_DIR}/software/angsd_dev/angsd"

echo -e "\nprocessing ${chunk}"

if [ -f "${outdir}/${string}.mafs.gz" ]; then
  echo "  Output already exists for ${chunk}, skipping."
  exit 0
fi

# mkdir is atomic, so only one job can take the lock
LOCKDIR="${outdir}/${string}.lock"
if ! mkdir "$LOCKDIR" 2>/dev/null; then
  echo "  ${chunk} is already being processed by another job, skipping."
  exit 0
fi
trap 'rmdir "$LOCKDIR" 2>/dev/null' EXIT

cd ${BASE_DIR}/mapping/merged

$ANGSD -P 9 -bam $BAMS -ref $REF -rf ${indir}/${string} -out ${outdir}/${string} \
-doGlf 2 -doCounts 1 -doMaf 1 -doPost 1 -baq 1 -doMajorMinor 4 -GL 2 -doBcf 1 \
-minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -skipTriallelic 1 -only_proper_pairs 0 \
-setMinDepthInd 1 -minInd 300 -setMinDepth 4623 -setMaxDepth 23116 -setMaxDepthInd 120 -doGeno 1 \
-SNP_pval 1e-7 -doHWE 1

# SNP FILTERS: 
# No MAF filter -- rare variants are critical for unbiased SFS and Tajima's D
# -SNP_pval 1e-7: keep only sites that are confidently variable
# -doGeno 1: required for the BCF GT field; -doHWE 1: per-site HWE estimates (no HWE filter applied)
# Paralogs are removed afterwards with ngsParalog (launch_ngsParalog.sh, site_filter.R)

# COVERAGE FILTERS:
# Individual depth ranges from ~1x to ~30x, in up to perhaps 980 plants.
# -minInd 300          sites must have data (depth >= setMinDepthInd) in at least 300 of 842 individuals;
#                      a site with fewer than 300 individuals at 1x almost certainly has a problem
# -setMinDepthInd 1    discard an individual's data below this depth (only for -doCounts analyses);
#                      kept at 1x because these individuals are combined with other datasets
# -setMinDepth 4623
# -setMaxDepth 23116   0.5x and 2.5x the summed mean autosomal depth across individuals
# -setMaxDepthInd 120  discard an individual's data above this depth: 2x the highest-coverage individual (60x)
