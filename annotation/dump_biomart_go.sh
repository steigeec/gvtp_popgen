#!/bin/bash
# Pull gene_id -> GO term pairs from Ensembl Plants BioMart for a species
# without a Bioconductor OrgDb, for build_ortholog_go.py. IEA evidence is excluded.
#
# Usage: bash dump_biomart_go.sh <biomart_dataset> <display_name>
#   biomart_dataset: Ensembl Plants BioMart dataset name, e.g. hannuus_eg_gene
#   display_name:    used for the output filename, e.g. Helianthus_annuus
#
# Example:
#   bash dump_biomart_go.sh hannuus_eg_gene Helianthus_annuus

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

if [ $# -lt 2 ]; then
    echo "Usage: bash dump_biomart_go.sh <biomart_dataset> <display_name>" >&2
    exit 1
fi

DATASET="$1"
DISPLAY_NAME="$2"

OUTDIR=${BASE_DIR}/genome/annotation
RAW="${OUTDIR}/${DISPLAY_NAME}_biomart_go_raw.tsv"
OUT="${OUTDIR}/${DISPLAY_NAME}_geneid_go.tsv"

QUERY=$(cat <<XML
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="plants_mart" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="${DATASET}" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="go_id"/>
    <Attribute name="go_linkage_type"/>
  </Dataset>
</Query>
XML
)

echo "Querying BioMart dataset ${DATASET}..."
curl -sS --fail --max-time 600 -G "https://plants.ensembl.org/biomart/martservice" \
    --data-urlencode "query=${QUERY}" -o "$RAW"

n_raw=$(wc -l < "$RAW")
echo "Received ${n_raw} raw (gene, GO, evidence) rows"

echo -e "gene_id\tgo_term" > "$OUT"
awk -F'\t' 'BEGIN{OFS="\t"} $1!="" && $2!="" && $3!="IEA" {print $1, $2}' "$RAW" | sort -u >> "$OUT"

n_pairs=$(($(wc -l < "$OUT") - 1))
n_genes=$(tail -n +2 "$OUT" | cut -f1 | sort -u | wc -l)
echo "Wrote ${n_pairs} unique gene-GO pairs across ${n_genes} genes (IEA evidence excluded) to ${OUT}"
