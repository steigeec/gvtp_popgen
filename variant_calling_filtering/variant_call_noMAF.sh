#!/bin/bash

chunk=$1
BAMS=$2
REF=$3
indir="/private/groups/meyerlab/esteiger/variant_calling/angsd/chromo_30_chunks/"
outdir="/private/groups/meyerlab/esteiger/variant_calling/angsd"
export TMPDIR="/tmp"
prefix="./"
string=${chunk#"$prefix"}
ANGSD="/private/groups/meyerlab/esteiger/software/angsd_dev/angsd"

echo -e "\nprocessing ${chunk}"

# Skip if output already exists
if [ -f "${outdir}/${string}.mafs.gz" ]; then
  echo "  Output already exists for ${chunk}, skipping."
  exit 0
fi

# Atomic lock using mkdir (only one process can create the directory)
LOCKDIR="${outdir}/${string}.lock"
if ! mkdir "$LOCKDIR" 2>/dev/null; then
  echo "  ${chunk} is already being processed by another job, skipping."
  exit 0
fi
trap 'rmdir "$LOCKDIR" 2>/dev/null' EXIT

cd /private/groups/meyerlab/esteiger/mapping/merged

$ANGSD -P 9 -bam $BAMS -ref $REF -rf ${indir}/${string} -out ${outdir}/${string} \
-doGlf 2 -doCounts 1 -doMaf 1 -doPost 1 -baq 1 -doMajorMinor 4 -GL 2 -doBcf 1 \
-minMapQ 20 -minQ 20 -remove_bads 1 -uniqueOnly 1 -skipTriallelic 1 -only_proper_pairs 0 \
-setMinDepthInd 1 -minInd 300 -setMinDepth 4623 -setMaxDepth 23116 -setMaxDepthInd 120 -doGeno 1 \
-SNP_pval 1e-7 -doHWE 1

# SNP FILTERS: 
# No MAF filter -- rare variants are critical for unbiased SFS and Tajima's D
# -doSnpStat 1: enables SNP artifact filters (requires ANGSD >= 0.941 for correct hetbias behaviour)
# -hetbias_pval 1e-5: removes sites with excess het in reads from het-called individuals (paralogs)
# -sb_pval 1e-4:   removes sites with strand bias in minor allele reads
# -edge_pval 1e-4: removes sites where minor allele reads cluster near read ends
# -mapQ_pval 1e-4: removes sites where minor allele reads have lower mapping quality
# -doHWE 1: required by doSnpStat; -doGeno 1: required for BCF GT field

# OUTPUTS:
# Beagle/BCF outputs included

# COVERAGE FILTERS:
# Should ranges from 1 to some 30x in up to perhaps 980 plants... 
# -minInd 300          <- There are 842 individuals total. Only keep sites with at least minIndDepth (default is 1) from at least [int] individuals. 
#	       If <300 individuals sequenced to 1x at that site, there's almost certainly something wrong with that site. 
# -setMinDepthInd 1    <- Discard individual if sequencing depth for an individual is below [int]. This filter is only applied to analysis which are based on counts of alleles i.e. analysis that uses -doCounts
#                        I go down to 1x here, because I'll be combining these individuals with my various datasets... 
# -setMinDepth 4623     
# -setMaxDepth 23116   <-between half and 2.5x the sum of the average autosomal sequencing depth of individuals 
# -setMaxDepthInd 120  <- Discard individual if sequencing depth for an individual is above [int]. 
#	       I use 2x my highest coverage individual, which is 60x 
