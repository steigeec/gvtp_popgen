#!/bin/bash
# Build a genome-wide SNP -> gene map (scaffold, position, gene_id) from the
# SnpEff-annotated VCF, for build_gene_gwas_scores.py. ANN parsing is in
# parse_snp_gene_ann.py; a SNP can map to more than one gene.
#
# Usage: bash build_snp_gene_map.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

VCF=${BASE_DIR}/selection/var_sites.ann.vcf
OUT=${BASE_DIR}/genome/annotation/snp_gene_map.tsv

echo "Extracting CHROM/POS/ANN from ${VCF}..."
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' "$VCF" \
    | python3 ${BASE_DIR}/gvtp_popgen/annotation/parse_snp_gene_ann.py \
    | sort -u -k1,1 -k2,2n -k3,3 > "$OUT"

echo "Wrote $(wc -l < "$OUT") (scaffold, position, gene_id) rows to ${OUT}"
echo "  Distinct SNPs: $(cut -f1,2 "$OUT" | sort -u | wc -l)"
echo "  Distinct genes: $(cut -f3 "$OUT" | sort -u | wc -l)"
