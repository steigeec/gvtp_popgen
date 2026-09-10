#!/bin/bash

#SBATCH --job-name=annot_40000_3
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=esteiger@ucsc.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=6:00:00
#SBATCH --output=annot_40000_3_%A_%a.log
#SBATCH --array=0-4

cd /private/groups/meyerlab/esteiger/selection/gemma/

VCF=/private/groups/meyerlab/esteiger/selection/var_sites.ann.vcf

# Putative inversion region — excluded from the GO enrichment candidate list
# (SNPs here are non-independent and would inflate any GO signal driven by
# the inversion itself rather than by selection at individual genes)
INV_SCAFFOLD="DEINVI_00006"
INV_START=43000000
INV_END=50000000

# "deficit" and "moisture" no longer have output directories (variables
# dropped from active investigation, per Launch_find_outlier_annotate.sh's
# own comment) -- kept this array in sync with that script's current
# active set instead. Harmless either way for BIO12/BCM_30y specifically:
# they have no outliers_for_annot_40000_3 file yet, so this script exits
# early for them regardless of what's in this list.
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
    # Extract gene list from filtered-peak annotations. Piped through
    # resolve_gene_ids.R instead of straight to sort -u/sed, since the raw
    # regex match isn't always already a gene ID (see that script's header).
    #
    # The awk step below extracts ONLY the Gene_ID subfield (pipe-field 5)
    # of each comma-separated ANN entry before the regex runs -- grepping
    # the raw ann.INFO column 5 directly (the old approach) scans every
    # pipe subfield at once, including Feature_ID, which can carry a raw
    # AUGUSTUS transcript stem like "g1478.t2" for a gene since renamed to
    # "IDmodified-gene-N". That stem matches the same regex and gets
    # counted as a spurious extra "g1478" gene call, silently colliding
    # with a real, unrelated gene elsewhere in the genome that kept id
    # "g1478" (confirmed: this is exactly how SCPL13 and UGT73E1 -- on
    # different scaffolds -- ended up labeled on the same BIO15 peak).
    awk -F'\t' 'NR>1 {n=split($5,e,","); for(i=1;i<=n;i++){split(e[i],p,"|"); print p[5]}}' \
        "${var_name}/outliers_40000_3_${var_name}.ann.INFO" | \
        grep -oE "MSTRG\.[0-9]+(\.p?[0-9]+)?|IDmodified-gene-[0-9]+|g[0-9]+|DEINVI_[0-9]+\.[0-9]+" | \
        sed "s/\(MSTRG\.[0-9]\+\)\.[0-9]\+.*/\1/" | \
        Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/resolve_gene_ids.R \
        > "${var_name}/gene_list_${var_name}_40000_3.txt"

    echo "Genes in 40000_3 set: $(wc -l < "${var_name}/gene_list_${var_name}_40000_3.txt")"

    # Exclude outlier SNPs that fall inside the putative inversion, then rebuild
    # the gene list as a subset of the full 40000_3 list. Genes whose only
    # supporting SNP(s) are inside the inversion drop out; genes with additional
    # non-inversion support are retained.
    awk -F'\t' -v OFS='\t' -v scaf="$INV_SCAFFOLD" -v s="$INV_START" -v e="$INV_END" \
        'NR==1 || !($1==scaf && $2>=s && $2<=e)' \
        "${var_name}/outliers_40000_3_${var_name}.ann.INFO" > "${var_name}/outliers_40000_3_${var_name}.ann.woInversion.INFO"

    # Same Gene_ID-subfield-only fix as above -- see that block's comment.
    awk -F'\t' 'NR>1 {n=split($5,e,","); for(i=1;i<=n;i++){split(e[i],p,"|"); print p[5]}}' \
        "${var_name}/outliers_40000_3_${var_name}.ann.woInversion.INFO" | \
        grep -oE "MSTRG\.[0-9]+(\.p?[0-9]+)?|IDmodified-gene-[0-9]+|g[0-9]+|DEINVI_[0-9]+\.[0-9]+" | \
        sed "s/\(MSTRG\.[0-9]\+\)\.[0-9]\+.*/\1/" | \
        Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/resolve_gene_ids.R \
        > "${var_name}/gene_list_${var_name}_40000_3_woInversion.txt"

    echo "Genes in 40000_3 set excluding inversion: $(wc -l < "${var_name}/gene_list_${var_name}_40000_3_woInversion.txt") (full set: $(wc -l < "${var_name}/gene_list_${var_name}_40000_3.txt"))"

    # Old GO enrichment (topGO/Fisher, no gene-length correction, includes inversion genes) — kept commented as a fallback
    # Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/GOenrich.R \
    #     "${var_name}/gene_list_${var_name}_40000_3.txt" "${var_name}"

    # GO enrichment for filtered peaks, controlling for gene length (goseq)
    Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/GOenrich_goseq.R \
        "${var_name}/gene_list_${var_name}_40000_3.txt" "${var_name}"

    # Only run it a second time on the woInversion list if the inversion
    # exclusion actually removed a gene — otherwise the two lists are
    # identical and a second run would just duplicate several minutes of
    # GO.db propagation for the same result (as happened with BIO18)
    if diff -q "${var_name}/gene_list_${var_name}_40000_3.txt" \
                "${var_name}/gene_list_${var_name}_40000_3_woInversion.txt" > /dev/null; then
        echo "No genes in the inversion for ${var_name} — skipping redundant woInversion GO run."
    else
        Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/GOenrich_goseq.R \
            "${var_name}/gene_list_${var_name}_40000_3_woInversion.txt" "${var_name}"
    fi

    # Tajima's D combined comparison plots for filtered peaks
    Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/selection/TajD/plot_TajD_extreme_combined.R \
        "40000_3"
else
    echo "No annotated variants in 40000_3 set for ${var_name} — skipping."
fi

echo "Done: ${var_name} 40000_3"
