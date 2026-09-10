#!/bin/bash
# Population-level allele frequencies for the AEM+RDA directional wind test
# (wind_aem_rda.R).
#
# ============================================================================
# PART 1 (PRIMARY): genotype-likelihood-based per-population MAF via ANGSD,
# reading GLs directly from the existing ANGSD-native BCF
# (noMAF_preParalog.bcf: angsd ... -GL 2 -doBcf 1, all 842 individuals, incl.
# down to 0.6x coverage) -- no BAM rescan needed, ANGSD's -vcf-PL reader just
# runs the -doMaf EM over the already-computed PL values.
#
# Site filtering: paralog-filtered (postParalog_auto.sites, 36.77M sites),
# but the existing pooled MAF>0.05 filter (postParalog_MAF_auto.sites,
# 7.94M sites) is deliberately NOT reapplied here. That pooled threshold is
# computed across all 842 individuals, of which Strauss is 738 (88%) -- so
# it is much closer to "common in Strauss" than a population-neutral filter,
# and would tend to exclude exactly the population-private/locally-elevated
# alleles (e.g. in Dangermond, n=4) that this analysis is testing for. The
# downstream drop_monomorphic() step in wind_aem_rda.R already filters on
# real variance *among the populations being compared*, which is the
# statistically appropriate target here, not pooled frequency across a
# hugely unbalanced pooled sample.
#
# doMajorMinor 4 (reference allele = major) matches -anc/-ref = same
# reference genome already used in pop_SAF.sh, keeping allele polarization
# consistent across all populations without needing a separate shared
# major/minor allele list.
# ============================================================================

set -euo pipefail

ANGSD_BIN=/private/groups/meyerlab/esteiger/software/angsd_dev/angsd
REF=/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked
ANGSD_DIR=/private/groups/meyerlab/esteiger/variant_calling/angsd
SITES="${ANGSD_DIR}/postParalog_auto.sites"     # paralog-filtered, NO pooled-MAF cut
FULL_BCF="${ANGSD_DIR}/noMAF_preParalog.bcf"    # ANGSD-native, all 842 individuals
SUBSET_BCF="${ANGSD_DIR}/postParalog_auto.bcf"  # NEW shared resource (written once, reused)

FIELD_DATA=/private/groups/meyerlab/esteiger/field_data.txt
GL_OUTDIR=/private/groups/meyerlab/esteiger/wind/allele_freq
mkdir -p "${GL_OUTDIR}"

# ---- Step 1: paralog-filtered, no-pooled-MAF-cut BCF (all 842 individuals) ----
if [ ! -f "${SUBSET_BCF}" ]; then
  echo "=== Building ${SUBSET_BCF} (subset of noMAF_preParalog.bcf to postParalog_auto.sites) ==="
  bcftools view -T "${SITES}" -O b -o "${SUBSET_BCF}" "${FULL_BCF}"
  bcftools index "${SUBSET_BCF}"
else
  echo "${SUBSET_BCF} already exists, skipping."
fi

# ---- Step 2: per-population sample lists ----
# BCF sample names are "<SampleID>_rmdup.bam"; one list per Clusters value.
awk -F'\t' '
  NR==1 { for(i=1;i<=NF;i++){ if($i=="SampleID") sid=i; if($i=="Clusters") cl=i }; next }
  { print $sid"_rmdup.bam" > ("'"${GL_OUTDIR}"'/samples_" $cl ".txt") }
' "${FIELD_DATA}"

# ---- Step 3: per-population ANGSD -doMaf, reading GLs directly from the BCF ----
for SAMPLE_LIST in "${GL_OUTDIR}"/samples_*.txt; do
  POP=$(basename "${SAMPLE_LIST}" .txt | sed 's/^samples_//')
  N=$(wc -l < "${SAMPLE_LIST}")
  POP_BCF="${GL_OUTDIR}/${POP}.bcf"
  echo "=== ${POP}: n=${N} individuals ==="

  bcftools view -S "${SAMPLE_LIST}" -O b -o "${POP_BCF}" "${SUBSET_BCF}"

  "${ANGSD_BIN}" -vcf-PL "${POP_BCF}" -nInd "${N}" -fai "${REF}.fai" -ref "${REF}" \
    -doMajorMinor 4 -doMaf 1 \
    -out "${GL_OUTDIR}/${POP}_maf"

  rm -f "${POP_BCF}"
  echo "  -> ${GL_OUTDIR}/${POP}_maf.mafs.gz"
done

echo "Done -- per-population GL-based MAF files in ${GL_OUTDIR}"

# ============================================================================
# PART 2 (ORIGINAL / superseded for this analysis): PLINK hard-called-genotype
# allele frequencies from the GEMMA dataset. Kept here for comparison only --
# Part 1 above is what wind_aem_rda.R should use, since hard-calling biases
# low-coverage heterozygotes and GEMMA's 4x filter drops all but 1 of
# Dangermond's 4 individuals (see chat discussion).
# ============================================================================

PLINK=/private/groups/meyerlab/esteiger/software/plink
BFILE=/private/groups/meyerlab/esteiger/selection/gemma/postParalog_MAF_4x_GLimp_clean
# Output (3+ GB) goes to env/, not the git-tracked project directory --
# matches where all other raw/derived wind data lives.
OUTDIR=/private/groups/meyerlab/esteiger/wind/allele_freq
TMP_BIM="/private/groups/meyerlab/esteiger/gvtp_popgen/wind/tmp_renamed.bim"
mkdir -p "${OUTDIR}"

awk 'BEGIN{OFS="\t"}{$1="s"$1; print}' "${BFILE}.bim" > "${TMP_BIM}"

"${PLINK}" \
  --bed "${BFILE}.bed" \
  --bim "${TMP_BIM}" \
  --fam "${BFILE}.fam" \
  --freq --family --allow-extra-chr \
  --out "${OUTDIR}/pop_allele_freq"

rm -f "${TMP_BIM}"
echo "Done -> ${OUTDIR}/pop_allele_freq.frq.strat"
