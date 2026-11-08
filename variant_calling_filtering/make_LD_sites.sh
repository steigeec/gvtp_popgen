#!/usr/bin/env bash
# Combine prune_graph .pos files into the LD-pruned sites file, then subset the
# beagle files to those sites: >4x individuals (PCAngsd PCs for GEMMA), all
# individuals, Year 1, and each locality.

BASE_DIR=/private/groups/meyerlab/esteiger

GITHUB="${BASE_DIR}/gvtp_popgen/variant_calling_filtering"
POS_DIR="${BASE_DIR}/variant_calling/ngsLD/pos"
SITES="${BASE_DIR}/variant_calling/ngsLD/angsd_paralog_LD.sites"
ANGSD_DIR="${BASE_DIR}/variant_calling/angsd"

# ---------- 1. LD-pruned sites ----------
n_pos=$(find "$POS_DIR" -name "DEINVI*_prune.pos" | wc -l)
echo "Found $n_pos .pos files"

cat "$POS_DIR"/DEINVI*_prune.pos \
    | sed 's/:/\t/' \
    | sort -k1,1 -k2,2n \
    | uniq \
    > "$SITES"
echo "Written: $SITES"
echo "Sites retained: $(wc -l < "$SITES")"

# ---------- 2. LD-pruned beagle files ----------

echo ""
echo "=== Create LD-filtered beagle subset to ONLY >4x individuals ==="
echo "=== To compute PCs for env covariates in GEMMA AFTER LD-filtering ... ==="

awk 'NR==FNR{keys[$1"_"$2]=1; next} FNR==1{print; next} ($1 in keys)' \
  "${SITES}" \
  <(zcat "${ANGSD_DIR}/postParalog_4x.beagle.gz") \
  | gzip > "${ANGSD_DIR}/postParalog_4x_LD.beagle.gz"&

# all individuals

awk 'NR==FNR{keys[$1"_"$2]=1; next} FNR==1{print; next} ($1 in keys)' \
  "${SITES}" \
  <(zcat "${ANGSD_DIR}/postParalog_MAF_auto.beagle.gz") \
  | gzip > "${ANGSD_DIR}/postParalog_LD.beagle.gz"

wait
echo "  -> $(zcat "${ANGSD_DIR}/postParalog_4x_LD.beagle.gz" | wc -l) lines (including header) in postParalog_4x_LD.beagle.gz"
echo "  -> $(zcat "${ANGSD_DIR}/postParalog_LD.beagle.gz" | wc -l) lines (including header) in postParalog_LD.beagle.gz"

# Year 1 and each locality
BEAGLE1="${ANGSD_DIR}/postParalog_LD.beagle.gz"
INDS="${BASE_DIR}/mapping/merged/bamlist_year1.txt"
OUT="angsd_paralog_LD_year1.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

INDS="${BASE_DIR}/mapping/merged/bamlist_strauss.txt"
OUT="angsd_paralog_LD_strauss.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

INDS="${BASE_DIR}/mapping/merged/bamlist_CDFW.txt"
OUT="angsd_paralog_LD_CDFW.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

INDS="${BASE_DIR}/mapping/merged/bamlist_Point_Conception.txt"
OUT="angsd_paralog_LD_Point_Conception.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

INDS="${BASE_DIR}/mapping/merged/bamlist_Dangermond.txt"
OUT="angsd_paralog_LD_Dangermond.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

INDS="${BASE_DIR}/mapping/merged/bamlist_VAFB_Lion_Head.txt"
OUT="angsd_paralog_LD_VAFB_Lion_Heads.beagle"
${GITHUB}/beagle_filter_individuals.sh ${INDS} ${BEAGLE1} ${OUT} &

wait

gzip angsd_paralog_LD_*.beagle

# Optional: LD-pruned imputed VCF (FEEMS input, migration/FEEMS.sh)
# bcftools view -R "${SITES}" \
#   "${ANGSD_DIR}/postParalog_MAF_auto_GLimp.vcf.gz" \
#   -O v -o "BASE_DIR/variant_calling/ngsLD/angsd_paralog_LD_GLimp.vcf"
