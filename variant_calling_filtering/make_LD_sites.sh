#!/usr/bin/env bash
# Concatenates all prune_graph .pos files, converts scaffold:pos to tab-separated,
# and writes a sorted, deduplicated sites file.
# Then creates the LD-filtered beagle file (>4x coverage sites) for PCAngsd
# computation of principal components used as covariates in GEMMA.

GITHUB="/private/groups/meyerlab/esteiger/gvtp_popgen/variant_calling_filtering"
POS_DIR="/private/groups/meyerlab/esteiger/variant_calling/ngsLD/pos"
SITES="/private/groups/meyerlab/esteiger/variant_calling/ngsLD/angsd_paralog_LD.sites"
ANGSD_DIR="/private/groups/meyerlab/esteiger/variant_calling/angsd"

# Step 1: build sites file
n_pos=$(find "$POS_DIR" -name "DEINVI*_prune.pos" | wc -l)
echo "Found $n_pos .pos files"

cat "$POS_DIR"/DEINVI*_prune.pos \
    | sed 's/:/\t/' \
    | sort -k1,1 -k2,2n \
    | uniq \
    > "$SITES"
echo "Written: $SITES"
echo "Sites retained: $(wc -l < "$SITES")"

############################################

echo ""
echo "=== Create LD-filtered beagle subset to ONLY >4x individuals ==="
echo "=== To compute PCs for env covariates in GEMMA AFTER LD-filtering ... ==="

awk 'NR==FNR{keys[$1"_"$2]=1; next} FNR==1{print; next} ($1 in keys)' \
  "${SITES}" \
  <(zcat "${ANGSD_DIR}/postParalog_4x.beagle.gz") \
  | gzip > "${ANGSD_DIR}/postParalog_4x_LD.beagle.gz"&

####################################
# Make LD-filtered beagle of ALL samples

awk 'NR==FNR{keys[$1"_"$2]=1; next} FNR==1{print; next} ($1 in keys)' \
  "${SITES}" \
  <(zcat "${ANGSD_DIR}/postParalog_MAF_auto.beagle.gz") \
  | gzip > postParalog_LD.beagle.gz"

wait
echo "  -> $(zcat postParalog_4x_LD.beagle.gz | wc -l) lines (including header) in postParalog_4x_LD.beagle.gz"
echo "  -> $(zcat postParalog_LD.beagle.gz | wc -l) lines (including header) in postParalog_LD.beagle.gz"

####################################
# Subsets for
# 	Year 1
#	each locality
BEAGLE1="${ANGSD_DIR}/postParalog_LD.beagle.gz"
INDS="/private/groups/meyerlab/esteiger/mapping/merged/bamlist_year1.txt"
OUT="angsd_paralog_LD_year1.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

INDS="/private/groups/meyerlab/esteiger/mapping/merged/bamlist_strauss.txt"
OUT="angsd_paralog_LD_strauss.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

INDS="/private/groups/meyerlab/esteiger/mapping/merged/bamlist_CDFW.txt"
OUT="angsd_paralog_LD_CDFW.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

INDS="/private/groups/meyerlab/esteiger/mapping/merged/bamlist_Point_Conception.txt"
OUT="angsd_paralog_LD_Point_Conception.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

INDS="/private/groups/meyerlab/esteiger/mapping/merged/bamlist_Dangermond.txt"
OUT="angsd_paralog_LD_Dangermond.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

INDS="/private/groups/meyerlab/esteiger/mapping/merged/bamlist_VAFB_Lion_Head.txt"
OUT="angsd_paralog_LD_VAFB_Lion_Heads.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

wait

gzip angsd_paralog_LD_*.beagle

####################################
# Creates the LD-filtered VCF (used e.g. as FEEMS input, see
# migration/FEEMS/FEEMS.sh) by subsetting postParalog_MAF_auto_GLimp.vcf to
# the LD-pruned site set above, via the tabix-indexed .vcf.gz for speed:
#
# bcftools view -R "${SITES}" \
#   "${ANGSD_DIR}/postParalog_MAF_auto_GLimp.vcf.gz" \
#   -O v -o "/private/groups/meyerlab/esteiger/variant_calling/ngsLD/angsd_paralog_LD_GLimp.vcf"
