#!/bin/bash
# ANGSD allele frequencies for one genome chunk, keeping invariant sites (input to Tajima's D).
# Usage: angsd_w_invar.sh <chunk> <bamlist> <reference.fasta>  (called by launch_angsd_w_invar.sh)

BASE_DIR=/private/groups/meyerlab/esteiger

chunk=$1
BAMS=$2
REF=$3
indir="${BASE_DIR}/variant_calling/angsd/chromo_30_chunks/"
outdir="${BASE_DIR}/selection/TajD"
export TMPDIR="/tmp"
prefix="./"
string=${chunk#"$prefix"}
ANGSD="${BASE_DIR}/software/angsd/angsd"

echo -e "\nprocessing ${chunk}"

# Rough coverage filters are applied here to throw out problematic (e.g. paralogous) sites early.
# Coverage: same filters as for variable sites.
# To keep all sites (variable and invariant): no min MAF filter, no SNP p-value filter, no HWE p-value filter.
# Other filters:
# -only_proper_pairs set to 0; I don't trust it anymore because of issues with polar bears.

cd ${BASE_DIR}/mapping/merged

$ANGSD -P 9 -bam $BAMS -ref $REF -rf ${indir}/${string} -out ${outdir}/${string} \
-doCounts 1 -doMaf 1 -doMajorMinor 4 -GL 2 -baq 1 \
-minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -skipTriallelic 1 -only_proper_pairs 0 \
-setMinDepthInd 1 -minInd 300 -setMinDepth 4623 -setMaxDepth 23116 -setMaxDepthInd 120

