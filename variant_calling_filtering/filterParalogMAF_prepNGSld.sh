#!/bin/bash
# After ngsParalog: build autosomal site lists (with and without MAF filter),
# the paralog/MAF-filtered beagle and BCF, per-scaffold ngsLD inputs for >4x
# Strauss individuals (then launches ngsLD_launch.sh), launch BEAGLE imputation,
# and a chloroplast site list.

BASE_DIR=/private/groups/meyerlab/esteiger

ANGSD_DIR="${BASE_DIR}/variant_calling/angsd"
PARA_DIR="${BASE_DIR}/variant_calling/ngsParalog"
ANGSD="${BASE_DIR}/software/angsd_dev/angsd"
MITO="${BASE_DIR}/genome/mito.txt"
CHLORO="${BASE_DIR}/genome/chloroplast.txt"
BAMLIST="${BASE_DIR}/mapping/merged/bamlist.txt"
BEAGLE="${ANGSD_DIR}/noMAF_preParalog.beagle.gz"
SCRIPTS_DIR="${BASE_DIR}/gvtp_popgen/variant_calling_filtering"

# ---------------------------------------------------------------------------
echo "=== Generate autosomal sites (exclude chloroplast & mitochondrial scaffolds) ==="
awk 'NR==FNR{exclude[$1]=1; next} !($1 in exclude) && $1 !~ /[0-9][A-Z]/' \
    <(cat "$CHLORO" "$MITO") \
    "${ANGSD_DIR}/postParalog_MAF.sites" \
    | sort -k1,1 -k2,2n \
    > "${ANGSD_DIR}/postParalog_MAF_auto.sites"
echo "  -> $(wc -l < ${ANGSD_DIR}/postParalog_MAF_auto.sites) autosomal sites in postParalog_MAF_auto.sites"

# ---------------------------------------------------------------------------
echo "=== Generate autosomal sites (no MAF filter; exclude chloroplast & mitochondrial scaffolds) ==="
# the regex drops malformed concatenated lines (e.g. DEINVI_00007DEINVI_00007);
# ANGSD sites index and tabix need position-sorted input
awk 'NR==FNR{exclude[$1]=1; next} !($1 in exclude) && $1 !~ /[0-9][A-Z]/' \
    <(cat "$CHLORO" "$MITO") \
    "${ANGSD_DIR}/postParalog.sites" \
    | sort -k1,1 -k2,2n \
    > "${ANGSD_DIR}/postParalog_auto.sites"
echo "  -> $(wc -l < ${ANGSD_DIR}/postParalog_auto.sites) autosomal sites in postParalog_auto.sites"
"$ANGSD" sites index "${ANGSD_DIR}/postParalog_auto.sites"
echo "  -> ANGSD index done"
cut -f1 "${ANGSD_DIR}/postParalog_auto.sites" | sort -u > "${ANGSD_DIR}/postParalog_auto.chrom"
echo "  -> $(wc -l < ${ANGSD_DIR}/postParalog_auto.chrom) chromosomes in postParalog_auto.chrom"

echo "=== Create beforeMAF autosomal VCF (for ngsF) ==="
bcftools view "${ANGSD_DIR}/noMAF_preParalog.bcf" \
  -T "${ANGSD_DIR}/postParalog_auto.sites" \
  -O v -o "${ANGSD_DIR}/postNGSparalog_beforeMAF_auto.vcf"
echo "  -> $(grep -cv '^#' ${ANGSD_DIR}/postNGSparalog_beforeMAF_auto.vcf) sites in postNGSparalog_beforeMAF_auto.vcf"

# ---------------------------------------------------------------------------
echo "=== Rename beagle header columns to sample names ==="
header=$(zcat "$BEAGLE" | head -n 1)
prefix=$(echo "$header" | cut -f1-3)
new_columns=$(awk '{for(i=1; i<=3; i++) printf "%s\t", $1}' "$BAMLIST" | sed 's/\t$//')
new_header=$(echo -e "$prefix\t$new_columns")
{
    echo -e "$new_header"
    zcat "$BEAGLE" | tail -n +2
} | gzip > "${ANGSD_DIR}/noMAF_preParalog_renamed.beagle.gz"
mv "${ANGSD_DIR}/noMAF_preParalog_renamed.beagle.gz" "$BEAGLE"
echo "  Done."

# ---------------------------------------------------------------------------
echo ""
echo "=== Create MAF & paralog-filtered beagle  ==="
echo "=== For various downstream processes  ==="

awk 'NR==FNR{keys[$1"_"$2]=1; next} FNR==1{print; next} ($1 in keys)' \
  "${ANGSD_DIR}/postParalog_MAF_auto.sites" \
  <(zcat "${ANGSD_DIR}/noMAF_preParalog_renamed.beagle.gz") \
  | gzip > "${ANGSD_DIR}/postParalog_MAF_auto.beagle.gz"
echo "  -> $(zcat ${ANGSD_DIR}/postParalog_MAF_auto.beagle.gz | wc -l) lines (including header) in postParalog_MAF_auto.beagle.gz"

# ---------------------------------------------------------------------------
echo ""
echo "=== Create beagle subset to ONLY Strauss individuals >4x ==="
echo "=== For LD computation and filtering... ==="

INDS_4X="${BASE_DIR}/mapping/merged/bamlist_4x_strauss.txt"
BEAGLE_FILTER="${SCRIPTS_DIR}/./beagle_filter_individuals.sh"
zcat "${ANGSD_DIR}/postParalog_MAF_auto.beagle.gz" > "${ANGSD_DIR}/postParalog_MAF_tmp.beagle" &
wait
"$BEAGLE_FILTER" "$INDS_4X" "${ANGSD_DIR}/postParalog_MAF_tmp.beagle" "${ANGSD_DIR}/postParalog_strauss_4x.beagle"
echo "  -> $(wc -l < ${ANGSD_DIR}/postParalog_strauss_4x.beagle) lines (including header) in postParalog_strauss_4x.beagle"
# all >4x individuals, for the PCs used as GEMMA covariates
INDS_4X="${BASE_DIR}/mapping/merged/bamlist_4x.txt"
"$BEAGLE_FILTER" "$INDS_4X" "${ANGSD_DIR}/postParalog_MAF_tmp.beagle" "${ANGSD_DIR}/postParalog_4x.beagle"
rm "${ANGSD_DIR}/postParalog_MAF_tmp.beagle"

NGSLD_DIR="${BASE_DIR}/variant_calling/ngsLD"
mkdir -p "$NGSLD_DIR"

echo "  Creating per-chromosome genotype likelihood TSVs for ngsLD..."
cat "${ANGSD_DIR}/postParalog_strauss_4x.beagle" | awk -v outdir="$NGSLD_DIR" '
NR > 1 {
    split($1, arr, "_");
    chrom = arr[1] "_" arr[2];
    $1=$2=$3="";
    print substr($0, 4) >> outdir "/" chrom ".tsv";
    close(outdir "/" chrom ".tsv");
}'
echo "  -> $(ls ${NGSLD_DIR}/*.tsv | wc -l) chromosome TSV files created"

echo "  Gzipping TSVs and launching ngsLD..."
cd "$NGSLD_DIR"
find . -name "*.tsv" | parallel -j 100 gzip
sbatch "${SCRIPTS_DIR}/ngsLD_launch.sh"

# ---------------------------------------------------------------------------
echo ""
echo "=== imputed VCFs ==="
echo "=== for various processes things ==="

bcftools view "${ANGSD_DIR}/noMAF_preParalog.bcf" \
  -T "${ANGSD_DIR}/postParalog_MAF_auto.sites" \
  -O b -o "${ANGSD_DIR}/postParalog_MAF_auto.bcf"
bcftools index "${ANGSD_DIR}/postParalog_MAF_auto.bcf"
echo "  -> $(bcftools view -H ${ANGSD_DIR}/postParalog_MAF_auto.bcf | wc -l) sites in postParalog_MAF_auto.bcf"

cd "$ANGSD_DIR"
sbatch ${SCRIPTS_DIR}/imputation.sh "${ANGSD_DIR}/postParalog_MAF_auto.bcf"

# >4x subset of the imputed VCF
FOUR="${BASE_DIR}/mapping/merged/bamlist_4x.txt"
VCF1="postParalog_MAF_auto_GLimp.vcf"
tabix -p vcf ${BASE_DIR}/variant_calling/angsd/postParalog_MAF_auto_GLimp.vcf.gz
bcftools view -S ${FOUR} -O v -o postParalog_MAF_4x_GLimp.vcf ${VCF1}
# remove invariable sites
bcftools view -c 1:minor postParalog_MAF_4x_GLimp.vcf -Ov -o postParalog_MAF_4x.variable.vcf
mv postParalog_MAF_4x.variable.vcf postParalog_MAF_4x_GLimp.vcf

# ---------------------------------------------------------------------------
echo ""
echo "=== Create chloroplast site list (paralog-filtered, no MAF filter) ==="
CHLORO_OUTDIR="${BASE_DIR}/structure/Fst/chloroplast"
mkdir -p "$CHLORO_OUTDIR"
grep -Ff "$CHLORO" "${PARA_DIR}/postParalog.sites" \
  > "${CHLORO_OUTDIR}/chloroplast.sites"
echo "  -> $(wc -l < ${CHLORO_OUTDIR}/chloroplast.sites) chloroplast sites"

# ---------------------------------------------------------------------------

