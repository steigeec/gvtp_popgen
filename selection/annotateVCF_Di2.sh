#!/bin/bash

#SBATCH --job-name=annotateVCF_Di2
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=esteiger@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time=168:00:00              # long partition max
#SBATCH --output=/private/groups/meyerlab/esteiger/gvtp_popgen/selection/annotateVCF_Di2_%j.log

# Re-annotates the outlier VCF against the corrected Di.2.0 SnpEff database
# (deduplicated, internal-stops-removed gene models — see
# /private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.gff).
# Only "ann" mode is run — "closest" mode (var_sites.closest.vcf) isn't used
# anywhere downstream, so it's skipped this time.
#
# Writes directly to var_sites.ann.vcf — this is now the canonical output;
# the pre-fix (Di.1.0) annotation is preserved at var_sites.ann.orig.vcf.
#
# Step 1 (below) strips the 842-sample genotype/FORMAT columns before running
# SnpEff: `ann` mode only uses CHROM/POS/REF/ALT/INFO to predict effects, it
# never reads genotypes, and every downstream step in this project only ever
# pulls the ANN INFO field back out (vcftools --get-INFO ANN) — the ~50-60KB
# of GT:DP:GL:PL:GP per site across 842 samples is unused weight that's the
# main reason this file is 1.8TB and SnpEff would otherwise take days.
# Skipped automatically if the stripped file already exists, so re-running
# this script after a failure downstream doesn't repeat the expensive part.

SNPEFF="/private/groups/meyerlab/esteiger/software/snpEff/snpEff.jar"
CONF="/private/groups/meyerlab/esteiger/software/snpEff/snpEff.config"

VCF=/private/groups/meyerlab/esteiger/variant_calling/angsd/postNGSparalog_beforeMAF_auto.vcf
SITES_ONLY=/private/groups/meyerlab/esteiger/selection/postNGSparalog_beforeMAF_auto.sites_only.vcf
OUT=/private/groups/meyerlab/esteiger/selection/var_sites.ann.vcf

if [ -s "$SITES_ONLY" ]; then
    echo "Sites-only VCF already exists at $SITES_ONLY — skipping strip step."
else
    echo "Stripping genotype columns from $VCF..."
    cut -f1-8 "$VCF" > "$SITES_ONLY"
    echo "Strip done. Line count (cut preserves line count exactly, so no need to re-read the 1.8TB original):"
    wc -l "$SITES_ONLY"
fi

STATS=/private/groups/meyerlab/esteiger/selection/snpEff_summary.html

# -stats pins the auto-generated summary (snpEff_summary.html) and gene-level
# stats (snpEff_genes.txt) next to $OUT — without it, SnpEff drops both into
# whatever directory this script happens to be sbatch'd from, which is easy
# to lose track of.
java -Xmx90g -jar "$SNPEFF" ann -c "$CONF" -v -stats "$STATS" -geneId Di.2.0 "$SITES_ONLY" > "$OUT"
