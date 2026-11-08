#!/bin/bash
# For each environmental variable (one per array task): GEMMA outliers (find_outliers.R),
# SnpEff annotation of outlier SNPs, LD-curated gene lists, orthologs and per-gene summary.
# Usage: sbatch Launch_find_outlier_annotate.sh  (STOP_BEFORE_GENE_LIST=1 stops after annotation)

#SBATCH --job-name=gemma_outliers
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=gemma_outlier_%A_%a.log
#SBATCH --array=0-4

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/selection/gemma/

VCF=${BASE_DIR}/selection/var_sites.ann.vcf
VCF_DS=${BASE_DIR}/variant_calling/angsd/postParalog_MAF_auto_GLimp_clean.vcf.gz
ORTHO=${BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.blastp.out
GO=${BASE_DIR}/genome/annotated_sequences_combined_geneid_go.tsv
GFF=${BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.gff
GENOME=${BASE_DIR}/genome/DEINVI_genome.fasta.masked

# Minimum r2 between an outlier SNP and a SNP inside a gene body to keep that gene when the
# SNP is annotated only as intergenic_region (needs ld_curation/ from launch_ld_scan.sh)
LD_R2_THRESHOLD=0.5

# Putative inversion region, excluded from the GO enrichment candidate list (SNPs here are non-independent)
INV_SCAFFOLD="DEINVI_00006"
INV_START=43000000
INV_END=50000000

# Column numbers in env.txt and variable names
pheno_nums=(3     4     5     6   11)
var_names=( BIO12 BIO15 BIO18 fog BCM_30y)

i=${SLURM_ARRAY_TASK_ID}
pheno_num="${pheno_nums[$i]}"
var_name="${var_names[$i]}"

echo "=========================================="
echo "Processing pheno${pheno_num}: ${var_name}"
echo "=========================================="

Rscript ${BASE_DIR}/gvtp_popgen/selection/gemma/find_outliers.R "$pheno_num" "$var_name"

outlier_file="${var_name}/unfiltered_outliers_for_annot_${var_name}.txt"

if [ ! -f "$outlier_file" ]; then
    echo "No outlier file found for ${var_name} — no outliers, done."
    exit 0
fi

if [ ! -s "$outlier_file" ]; then
    echo "No outliers for ${var_name} — done."
    exit 0
fi

# Optional: extract imputed dosages (DS) for outlier SNPs (dosage_<var>.txt)
# if [ -f "${var_name}/dosage_${var_name}.txt" ]; then
#     echo "Dosage file already exists for ${var_name} — skipping extraction."
# else
#     echo "Extracting dosage values for ${var_name}..."
#     ds_outlier_file="${var_name}/outliers_for_annot_40000_3_${var_name}.txt"
#     if [ ! -f "$ds_outlier_file" ]; then
#         ds_outlier_file="$outlier_file"
#     fi
#
#     awk '{print $2}' "$ds_outlier_file" | sort -u > "${var_name}/tmp_snps.txt"
#     awk 'NR==FNR {s[$1]; next} $2 in s {print $3"\t"$4"\t"$2}' \
#         "${var_name}/tmp_snps.txt" scaffold_snp_map.txt > "${var_name}/tmp_snp_regions.txt"
#     awk '{print $1"\t"$2}' "${var_name}/tmp_snp_regions.txt" > "${var_name}/tmp_regions.txt"
#
#     printf "snp\t" > "${var_name}/dosage_${var_name}.txt"
#     bcftools query -l "$VCF_DS" | sed 's/_rmdup\.bam//' | paste -sd'\t' >> "${var_name}/dosage_${var_name}.txt"
#     bcftools query -R "${var_name}/tmp_regions.txt" -f '%CHROM\t%POS[\t%DS]\n' "$VCF_DS" \
#         | awk 'NR==FNR {map[$1"\t"$2]=$3; next}
#                {key=$1"\t"$2; if(key in map) {printf "%s", map[key]; for(i=3;i<=NF;i++) printf "\t%s", $i; printf "\n"}}' \
#             "${var_name}/tmp_snp_regions.txt" - >> "${var_name}/dosage_${var_name}.txt"
#     rm "${var_name}/tmp_snps.txt" "${var_name}/tmp_snp_regions.txt" "${var_name}/tmp_regions.txt"
#     echo "Dosage file: $(wc -l < "${var_name}/dosage_${var_name}.txt") lines"
# fi

echo "Annotating ${var_name}..."

# Map internal scaffold+SNP IDs to VCF chromosome+position (unfiltered)
cut -f1,2 -d" " "$outlier_file" | \
awk -v OFS="\t" '
    NR==FNR {
        map[$1"\t"$2] = $3"\t"$4
        next
    }
    {
        key = $1"\t"$2
        if (key in map) print map[key]
    }' scaffold_snp_map.txt - > "${var_name}/snps_${var_name}.txt"

# Extract variant annotations from VCF (unfiltered)
vcftools --vcf "$VCF" \
         --out "${var_name}/outliers_${var_name}.ann" \
         --positions "${var_name}/snps_${var_name}.txt" \
         --get-INFO ANN

ann_lines=$(tail -n +2 "${var_name}/outliers_${var_name}.ann.INFO" | wc -l)
if [ "$ann_lines" -eq 0 ]; then
    echo "No annotated variants found for ${var_name} — skipping downstream annotation."
    exit 0
fi

# Annotate the 40000_3 (retained-peak) SNP set
filtered_outlier_file="${var_name}/outliers_for_annot_40000_3_${var_name}.txt"
have_40000_3_ann=0
if [ -f "$filtered_outlier_file" ] && [ -s "$filtered_outlier_file" ]; then
    echo "Annotating 40000_3 peak set for ${var_name}..."

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

    vcftools --vcf "$VCF" \
             --out "${var_name}/outliers_40000_3_${var_name}.ann" \
             --positions "${var_name}/snps_40000_3_${var_name}.txt" \
             --get-INFO ANN

    ann40_lines=$(tail -n +2 "${var_name}/outliers_40000_3_${var_name}.ann.INFO" | wc -l)
    if [ "$ann40_lines" -gt 0 ]; then
        have_40000_3_ann=1
    else
        echo "No annotated variants in 40000_3 set for ${var_name}."
    fi
else
    echo "No 40000_3 outlier file for ${var_name} — skipping filtered-peak annotation."
fi

# For a new SNP set: stop here, run find_intergenic_candidate_genes.py on the .ann.INFO files
# and launch_ld_scan.sh, then re-run without the flag.
if [ "${STOP_BEFORE_GENE_LIST:-0}" = "1" ]; then
    echo "STOP_BEFORE_GENE_LIST=1 — stopping after annotation for ${var_name}."
    exit 0
fi

# Gene list (unfiltered): intergenic-only calls kept only with LD support (curate_intergenic_by_ld.py);
# resolve_gene_ids.R maps transcript IDs to gene IDs.
python3 ${BASE_DIR}/gvtp_popgen/selection/gemma/curate_intergenic_by_ld.py \
    "${var_name}/outliers_${var_name}.ann.INFO" \
    --ld-curation-dir ${BASE_DIR}/selection/gemma/ld_curation \
    --r2-threshold "$LD_R2_THRESHOLD" | \
    Rscript ${BASE_DIR}/gvtp_popgen/selection/gemma/resolve_gene_ids.R \
    > "${var_name}/gene_list_${var_name}_unfiltered.txt"

# Gene list (40000_3), used by compile_gene_summary.py to flag genes in retained peaks
if [ "$have_40000_3_ann" -eq 1 ]; then
    python3 ${BASE_DIR}/gvtp_popgen/selection/gemma/curate_intergenic_by_ld.py \
        "${var_name}/outliers_40000_3_${var_name}.ann.INFO" \
        --ld-curation-dir ${BASE_DIR}/selection/gemma/ld_curation \
        --r2-threshold "$LD_R2_THRESHOLD" | \
        Rscript ${BASE_DIR}/gvtp_popgen/selection/gemma/resolve_gene_ids.R \
        > "${var_name}/gene_list_${var_name}_40000_3.txt"
fi

# Orthologs (Asteraceae BLAST, then OrthoFinder Arabidopsis, then eggNOG; see build_gene_orthology.py)
python3 ${BASE_DIR}/gvtp_popgen/selection/gemma/build_gene_orthology.py \
    "${var_name}/gene_list_${var_name}_unfiltered.txt" --out "${var_name}/orthologs_${var_name}.txt"

# Optional: GO-term lookup for the gene summary (needs $GO from annotation/merge_go_sources.sh)
# python3 BASE_DIR/gvtp_popgen/selection/gemma/lookup_by_gene_id.py \
#     --by-gene "${var_name}/gene_list_${var_name}_unfiltered.txt" "$GO" > "${var_name}/go_${var_name}.txt"

# Per-gene summary (gene list, orthologs, GO terms if present, annotations, 40000_3 flag)
python3 ${BASE_DIR}/gvtp_popgen/selection/gemma/compile_gene_summary.py "${var_name}"

# Optional: 40000_3 peak plots and analyses
# Rscript BASE_DIR/gvtp_popgen/selection/gemma/annotated_manhattan.R "$pheno_num" "$var_name"
#
# if [ "$var_name" = "fog" ] || [ "$var_name" = "BIO15" ]; then
#     Rscript BASE_DIR/gvtp_popgen/selection/gemma/fog_deinvi6_figure.R "$var_name"
# fi
#
# Rscript BASE_DIR/gvtp_popgen/selection/gemma/germ_year_polygenic.R \
#     "$var_name" "$var_name" "${var_name}_germ_yr"

echo "Done: ${var_name}"
