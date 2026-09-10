#!/bin/bash

chunk=$1
BAMS=$2
REF=$3
indir="/private/groups/meyerlab/esteiger/variant_calling/angsd/chromo_30_chunks/"
outdir="/private/groups/meyerlab/esteiger/selection/TajD"
export TMPDIR="/tmp"
prefix="./"
string=${chunk#"$prefix"}
ANGSD="/private/groups/meyerlab/esteiger/software/angsd/angsd"

echo -e "\nprocessing ${chunk}"

# Here we do want to apply some rough coverage filters, which will help us early on throw out problematic (e.g. paralogous) sites
# COVERAGE... 
# 	same as for variable sites... 
# TO KEEP ALL VARIABLE
#	no min maf filter
# 	no SNP pvalue filter
#	no HWE pvalue filter

# OTHER FILTERS
# -only_proper_pairs I don't trust anymore, because of issues with polar bears, so set to 0

cd /private/groups/meyerlab/esteiger/mapping/merged

$ANGSD -P 9 -bam $BAMS -ref $REF -rf ${indir}/${string} -out ${outdir}/${string} \
-doCounts 1 -doMaf 1 -doMajorMinor 4 -GL 2 -baq 1 \
-minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -skipTriallelic 1 -only_proper_pairs 0 \
-setMinDepthInd 1 -minInd 300 -setMinDepth 4623 -setMaxDepth 23116 -setMaxDepthInd 120

