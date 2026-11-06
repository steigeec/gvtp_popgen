#!/bin/bash

#SBATCH --job-name=gemma_outliers
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=gemma_outlier_%A_%a.log
#SBATCH --array=0-4

cd /private/groups/meyerlab/esteiger/selection/gemma/

VCF=/private/groups/meyerlab/esteiger/selection/var_sites.ann.vcf
VCF_DS=/private/groups/meyerlab/esteiger/variant_calling/angsd/postParalog_MAF_auto_GLimp_clean.vcf.gz
ORTHO=/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.blastp.out
GO=/private/groups/meyerlab/esteiger/genome/annotated_sequences_combined_geneid_go.tsv
GFF=/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.gff
GENOME=/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked

# Minimum LD (r2) between an outlier SNP and a SNP inside a candidate gene's
# body required to keep that gene when the SNP's ONLY snpEff annotation is
# intergenic_region (unbounded nearest-flanking-gene distance). Requires the
# one-time genome-wide precompute in ld_curation/ — see
# find_intergenic_candidate_genes.py and launch_ld_scan.sh.
LD_R2_THRESHOLD=0.5

# Putative inversion region — excluded from the GO enrichment candidate list
# (SNPs here are non-independent and would inflate any GO signal driven by
# the inversion itself rather than by selection at individual genes)
INV_SCAFFOLD="DEINVI_00006"
INV_START=43000000
INV_END=50000000

# Restricted to the variables still under active investigation (their
# output directories are the only ones that survived cleanup of the
# original 16-variable sweep — the rest are no longer of interest).
pheno_nums=(3     4     5     6   11)
var_names=( BIO12 BIO15 BIO18 fog BCM_30y)

i=${SLURM_ARRAY_TASK_ID}
pheno_num="${pheno_nums[$i]}"
var_name="${var_names[$i]}"

echo "=========================================="
echo "Processing pheno${pheno_num}: ${var_name}"
echo "=========================================="

Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/find_outliers.R "$pheno_num" "$var_name"

# Run TajD combined plots once (array task 0 only — script reads all variables internally)
# if [ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]; then
#     Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/TajD/plot_TajD_extreme_combined.R
#     Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/TajD/plot_TajD_extreme_combined.R "40000_3"
# fi

outlier_file="${var_name}/unfiltered_outliers_for_annot_${var_name}.txt"

if [ ! -f "$outlier_file" ]; then
    echo "No outlier file found for ${var_name} — no outliers, done."
    exit 0
fi

if [ ! -s "$outlier_file" ]; then
    echo "No outliers for ${var_name} — done."
    exit 0
fi

# # Extract dosage (DS) values from imputed VCF for outlier SNPs
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

# Also annotate the 40000_3 (retained-peak) SNP set here, ahead of any
# gene-list creation — so both ann.INFO files are on disk before the
# LD-informed gene-list step runs once, instead of annotating 40000_3 late
# (as this script used to) and needing a second gene-list pass.
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

# Stop here when regenerating annotation for a fresh SNP set whose scaffolds
# the LD scan (ld_curation/) hasn't covered yet. Run
# find_intergenic_candidate_genes.py against the fresh .ann.INFO files, then
# launch_ld_scan.sh, then re-run this script without the flag so gene-list
# creation below runs exactly once with full LD coverage instead of twice.
if [ "${STOP_BEFORE_GENE_LIST:-0}" = "1" ]; then
    echo "STOP_BEFORE_GENE_LIST=1 — stopping after annotation for ${var_name}."
    exit 0
fi

# Extract gene list (unfiltered). Genic and bounded upstream/downstream_gene_variant
# calls are kept unconditionally; a SNP whose ONLY annotation is
# intergenic_region (unbounded nearest-flanking-gene, no distance cap) only
# keeps that gene if outlier SNP and gene body share LD r2 >= LD_R2_THRESHOLD
# per the precomputed ld_curation/ scan — see curate_intergenic_by_ld.py's
# header. Falls back to dropping (not keeping) any intergenic-only call on a
# scaffold the LD scan hasn't covered yet — run find_intergenic_candidate_genes.py
# + launch_ld_scan.sh first (one-time, genome-wide) for full coverage.
# Piped through resolve_gene_ids.R instead of straight to sort -u/sed, since
# the raw regex match isn't always already a gene ID (some AUGUSTUS/
# transdecoder transcript ids are numbered independently from their parent
# gene id — see that script's header for the concrete example).
python3 /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/curate_intergenic_by_ld.py \
    "${var_name}/outliers_${var_name}.ann.INFO" \
    --ld-curation-dir /private/groups/meyerlab/esteiger/selection/gemma/ld_curation \
    --r2-threshold "$LD_R2_THRESHOLD" | \
    Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/resolve_gene_ids.R \
    > "${var_name}/gene_list_${var_name}_unfiltered.txt"

# Extract gene list (40000_3) too, same LD-informed curation — kept only so
# compile_gene_summary.py can flag which unfiltered genes also fall in a
# 40000_3 peak. The GO-enrichment/inversion machinery for this scope stays
# commented out below along with the unfiltered one (see that block).
if [ "$have_40000_3_ann" -eq 1 ]; then
    python3 /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/curate_intergenic_by_ld.py \
        "${var_name}/outliers_40000_3_${var_name}.ann.INFO" \
        --ld-curation-dir /private/groups/meyerlab/esteiger/selection/gemma/ld_curation \
        --r2-threshold "$LD_R2_THRESHOLD" | \
        Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/resolve_gene_ids.R \
        > "${var_name}/gene_list_${var_name}_40000_3.txt"
fi

# Raw BLASTP top-hit (any species) isn't orthology inference -- it produced
# spurious cross-kingdom "orthologs" during the gene_summary_all.tsv
# literature pass (e.g. a mammalian sperm-specific gene as the reported
# ortholog). build_gene_orthology.py instead prioritizes real orthology
# calls: a strict-evalue sunflower-specific BLAST hit, then OrthoFinder's
# Arabidopsis ortholog call (described via the local TAIR10 FASTA), then
# eggNOG's own description for the gene -- falling back to the general
# BLAST top hit (as BLAST_Gene_Names, not Orthologs) only when none of
# those found anything. See that script's header for the full priority
# chain and reasoning.
python3 /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/build_gene_orthology.py \
    "${var_name}/gene_list_${var_name}_unfiltered.txt" --out "${var_name}/orthologs_${var_name}.txt"

# ── GO-dependent steps: commented out while functional annotation ($GO) is
# being rebuilt. Covers GO-term harvest, the inversion/woInversion gene-list
# rebuild (it only exists to give GO enrichment a with/without-inversion
# comparison), and GO enrichment itself, for both the unfiltered and 40000_3
# scopes. Uncomment once the new $GO file is ready.
#
# python3 /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/lookup_by_gene_id.py \
#     --by-gene "${var_name}/gene_list_${var_name}_unfiltered.txt" "$GO" > "${var_name}/go_${var_name}.txt"
#
# awk -F'\t' -v OFS='\t' -v scaf="$INV_SCAFFOLD" -v s="$INV_START" -v e="$INV_END" \
#     'NR==1 || !($1==scaf && $2>=s && $2<=e)' \
#     "${var_name}/outliers_${var_name}.ann.INFO" > "${var_name}/outliers_${var_name}.ann.woInversion.INFO"
#
# python3 /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/curate_intergenic_by_ld.py \
#     "${var_name}/outliers_${var_name}.ann.woInversion.INFO" \
#     --ld-curation-dir /private/groups/meyerlab/esteiger/selection/gemma/ld_curation \
#     --r2-threshold "$LD_R2_THRESHOLD" | \
#     Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/resolve_gene_ids.R \
#     > "${var_name}/gene_list_${var_name}_unfiltered_woInversion.txt"
#
# Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/GOenrich_goseq.R \
#     "${var_name}/gene_list_${var_name}_unfiltered.txt" "${var_name}"
#
# Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/top_nominal_GO_hits.R \
#     "${var_name}/gene_list_${var_name}_unfiltered.txt" "${var_name}"
#
# if diff -q "${var_name}/gene_list_${var_name}_unfiltered.txt" \
#             "${var_name}/gene_list_${var_name}_unfiltered_woInversion.txt" > /dev/null; then
#     echo "No genes in the inversion for ${var_name} — skipping redundant woInversion GO run."
# else
#     Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/GOenrich_goseq.R \
#         "${var_name}/gene_list_${var_name}_unfiltered_woInversion.txt" "${var_name}"
#     Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/top_nominal_GO_hits.R \
#         "${var_name}/gene_list_${var_name}_unfiltered_woInversion.txt" "${var_name}"
# fi
#
# if [ "$have_40000_3_ann" -eq 1 ]; then
#     awk -F'\t' -v OFS='\t' -v scaf="$INV_SCAFFOLD" -v s="$INV_START" -v e="$INV_END" \
#         'NR==1 || !($1==scaf && $2>=s && $2<=e)' \
#         "${var_name}/outliers_40000_3_${var_name}.ann.INFO" > "${var_name}/outliers_40000_3_${var_name}.ann.woInversion.INFO"
#
#     python3 /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/curate_intergenic_by_ld.py \
#         "${var_name}/outliers_40000_3_${var_name}.ann.woInversion.INFO" \
#         --ld-curation-dir /private/groups/meyerlab/esteiger/selection/gemma/ld_curation \
#         --r2-threshold "$LD_R2_THRESHOLD" | \
#         Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/resolve_gene_ids.R \
#         > "${var_name}/gene_list_${var_name}_40000_3_woInversion.txt"
#
#     Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/GOenrich_goseq.R \
#         "${var_name}/gene_list_${var_name}_40000_3.txt" "${var_name}"
#
#     Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/top_nominal_GO_hits.R \
#         "${var_name}/gene_list_${var_name}_40000_3.txt" "${var_name}"
#
#     if diff -q "${var_name}/gene_list_${var_name}_40000_3.txt" \
#                 "${var_name}/gene_list_${var_name}_40000_3_woInversion.txt" > /dev/null; then
#         echo "No genes in the inversion for ${var_name} — skipping redundant woInversion GO run."
#     else
#         Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/GOenrich_goseq.R \
#             "${var_name}/gene_list_${var_name}_40000_3_woInversion.txt" "${var_name}"
#         Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/top_nominal_GO_hits.R \
#             "${var_name}/gene_list_${var_name}_40000_3_woInversion.txt" "${var_name}"
#     fi
# fi

# Compile per-gene summary from the unfiltered gene list, orthologs, GO
# terms (if present), and outlier annotations produced by the steps above.
# Also flags whether each gene falls in a 40000_3 peak, via gene_list_*_40000_3.txt.
python3 /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/compile_gene_summary.py "${var_name}"

# ── 40000_3-dependent plots/analyses: commented out along with the rest of
# the 40000_3 scope (see note above the 40000_3 gene-list extraction).
# Uncomment once the 40000_3 branch is back in active use.
#
# Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/annotated_manhattan.R "$pheno_num" "$var_name"
#
# if [ "$var_name" = "fog" ] || [ "$var_name" = "BIO15" ]; then
#     Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/fog_deinvi6_figure.R "$var_name"
# fi
#
# Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/germ_year_polygenic.R \
#     "$var_name" "$var_name" "${var_name}_germ_yr"

echo "Done: ${var_name}"
