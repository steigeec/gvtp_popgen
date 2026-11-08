#!/bin/bash
# Annotate variant sites with SnpEff against the Di.2.0 database
# (gene models: BASE_DIR/genome/remove_duplicates_rm_internalstops_longest_isoform.gff).
# Output: var_sites.ann.vcf and snpEff_summary.html in BASE_DIR/selection

#SBATCH --job-name=annotateVCF_Di2
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=168:00:00              # long partition max
#SBATCH --output=/private/groups/meyerlab/esteiger/gvtp_popgen/selection/annotateVCF_Di2_%j.log

BASE_DIR=/private/groups/meyerlab/esteiger

SNPEFF="${BASE_DIR}/software/snpEff/snpEff.jar"
CONF="${BASE_DIR}/software/snpEff/snpEff.config"

VCF=${BASE_DIR}/variant_calling/angsd/postNGSparalog_beforeMAF_auto.vcf
SITES_ONLY=${BASE_DIR}/selection/postNGSparalog_beforeMAF_auto.sites_only.vcf
OUT=${BASE_DIR}/selection/var_sites.ann.vcf

if [ -s "$SITES_ONLY" ]; then
    echo "Sites-only VCF already exists at $SITES_ONLY — skipping strip step."
else
    echo "Stripping genotype columns from $VCF..."
    cut -f1-8 "$VCF" > "$SITES_ONLY"
    echo "Strip done. Line count (cut preserves line count exactly, so no need to re-read the 1.8TB original):"
    wc -l "$SITES_ONLY"
fi

# Genotype columns are dropped: SnpEff ann only needs the first 8 columns, and downstream steps read only ANN.

STATS=${BASE_DIR}/selection/snpEff_summary.html

# -stats writes snpEff_summary.html and snpEff_genes.txt next to $OUT rather than the submission directory.
java -Xmx90g -jar "$SNPEFF" ann -c "$CONF" -v -stats "$STATS" -geneId Di.2.0 "$SITES_ONLY" > "$OUT"
