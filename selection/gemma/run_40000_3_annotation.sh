#!/bin/bash
# SnpEff annotations and gene lists (with and without the putative inversion) for the
# 40000_3 filtered outlier SNPs of each environmental variable, one per array task.

#SBATCH --job-name=annot_40000_3
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=6:00:00
#SBATCH --output=annot_40000_3_%A_%a.log
#SBATCH --array=0-4

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/selection/gemma/

VCF=${BASE_DIR}/selection/var_sites.ann.vcf

# Putative inversion region, excluded from the GO enrichment candidate list (SNPs here are non-independent)
INV_SCAFFOLD="DEINVI_00006"
INV_START=43000000
INV_END=50000000

# Same variables as Launch_find_outlier_annotate.sh; variables without an outlier file exit early.
var_names=(BIO12 BIO15 BIO18 fog BCM_30y)

i=${SLURM_ARRAY_TASK_ID}
var_name="${var_names[$i]}"

echo "=========================================="
echo "40000_3 annotation: ${var_name}"
echo "=========================================="

filtered_outlier_file="${var_name}/outliers_for_annot_40000_3_${var_name}.txt"

if [ ! -f "$filtered_outlier_file" ] || [ ! -s "$filtered_outlier_file" ]; then
    echo "No 40000_3 outlier file for ${var_name} — exiting."
    exit 0
fi

# Map filtered outlier SNPs to VCF chromosome+position
cut -f1,2 -d" " "$filtered_outlier_file" | \
awk -v OFS="\t" '
    NR==FNR {
        map[$1"\t"$2] = $3"\t"$4
        next
    }
    {
        key = $1"\t"$2
        if (key in map) print map[key]
    }' scaffold_snp_map.txt - > "${var_name}/snps_40000_3_${var_name}.txt"

echo "SNP positions extracted: $(wc -l < "${var_name}/snps_40000_3_${var_name}.txt")"

# Extract variant annotations from VCF for filtered peaks
vcftools --vcf "$VCF" \
         --out "${var_name}/outliers_40000_3_${var_name}.ann" \
         --positions "${var_name}/snps_40000_3_${var_name}.txt" \
         --get-INFO ANN

ann40_lines=$(tail -n +2 "${var_name}/outliers_40000_3_${var_name}.ann.INFO" | wc -l)
echo "Annotated variants: ${ann40_lines}"

if [ "$ann40_lines" -gt 0 ]; then
    # Gene list from the ANN Gene_ID subfield only (Feature_ID can carry old AUGUSTUS
    # "g<number>.t<N>" stems that match unrelated genes); resolve_gene_ids.R maps
    # transcript IDs to gene IDs.
    awk -F'\t' 'NR>1 {n=split($5,e,","); for(i=1;i<=n;i++){split(e[i],p,"|"); print p[5]}}' \
        "${var_name}/outliers_40000_3_${var_name}.ann.INFO" | \
        grep -oE "MSTRG\.[0-9]+(\.p?[0-9]+)?|IDmodified-gene-[0-9]+|g[0-9]+|DEINVI_[0-9]+\.[0-9]+" | \
        sed "s/\(MSTRG\.[0-9]\+\)\.[0-9]\+.*/\1/" | \
        Rscript ${BASE_DIR}/gvtp_popgen/selection/gemma/resolve_gene_ids.R \
        > "${var_name}/gene_list_${var_name}_40000_3.txt"

    echo "Genes in 40000_3 set: $(wc -l < "${var_name}/gene_list_${var_name}_40000_3.txt")"

    # Gene list again after excluding outlier SNPs inside the putative inversion
    awk -F'\t' -v OFS='\t' -v scaf="$INV_SCAFFOLD" -v s="$INV_START" -v e="$INV_END" \
        'NR==1 || !($1==scaf && $2>=s && $2<=e)' \
        "${var_name}/outliers_40000_3_${var_name}.ann.INFO" > "${var_name}/outliers_40000_3_${var_name}.ann.woInversion.INFO"

    awk -F'\t' 'NR>1 {n=split($5,e,","); for(i=1;i<=n;i++){split(e[i],p,"|"); print p[5]}}' \
        "${var_name}/outliers_40000_3_${var_name}.ann.woInversion.INFO" | \
        grep -oE "MSTRG\.[0-9]+(\.p?[0-9]+)?|IDmodified-gene-[0-9]+|g[0-9]+|DEINVI_[0-9]+\.[0-9]+" | \
        sed "s/\(MSTRG\.[0-9]\+\)\.[0-9]\+.*/\1/" | \
        Rscript ${BASE_DIR}/gvtp_popgen/selection/gemma/resolve_gene_ids.R \
        > "${var_name}/gene_list_${var_name}_40000_3_woInversion.txt"

    echo "Genes in 40000_3 set excluding inversion: $(wc -l < "${var_name}/gene_list_${var_name}_40000_3_woInversion.txt") (full set: $(wc -l < "${var_name}/gene_list_${var_name}_40000_3.txt"))"
else
    echo "No annotated variants in 40000_3 set for ${var_name} — skipping."
fi

echo "Done: ${var_name} 40000_3"
