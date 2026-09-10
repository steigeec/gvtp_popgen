#!/bin/bash
#SBATCH --job-name=per_cohort_Ne
#SBATCH --partition=long
#SBATCH --cpus-per-task=15
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/demography/ANGSD-Ne/logs/%x.out
#SBATCH --error=/private/groups/meyerlab/esteiger/demography/ANGSD-Ne/logs/%x.err

set -euo pipefail

# Per-cohort (per-sampling-year) Ne for Tranquillon Mtn (strauss), comparing three
# sampling years sampled from the same population. SAF files (strauss_year1/2/3) are
# built by pop_SAF.sh from the year-specific >4x-depth bamlists
# (bamlist_year{1,2,3}_4x_strauss.txt), written to SAFDIR below.
#
# This mirrors the population-level pipeline (pop_het.sh -> run_watterson_theta.sh ->
# summarize_watterson_theta.sh) but folds all three steps into one script since there
# are only 3 cohorts to process.

MISC=/private/groups/meyerlab/esteiger/software/angsd_dev/misc
RS=/private/groups/meyerlab/esteiger/software/angsd/misc/realSFS
SAFDIR=/private/groups/meyerlab/esteiger/SAF/population
OUTDIR=/private/groups/meyerlab/esteiger/demography/ANGSD-Ne/per-cohort
CORES=${SLURM_CPUS_PER_TASK:-15}
MU=1e-8   # per site per generation; Asteraceae rate (Strasburg & Rieseberg 2008,
          # Helianthus annuus/petiolaris; Senecio chrysanthemifolius, ~1e-8) -- matches
          # summarize_watterson_theta.sh so Ne estimates are directly comparable across years

COHORTS=(strauss_year1 strauss_year2 strauss_year3)

mkdir -p "$OUTDIR" "$OUTDIR/logs"
cd "$OUTDIR"

SUMMARY="${OUTDIR}/per_cohort_watterson_theta_Ne_summary.tsv"
echo -e "cohort\tsum_tW\tsum_nSites\ttheta_W_per_site\tNe" > "$SUMMARY"

for COHORT in "${COHORTS[@]}"; do
    SAF="${SAFDIR}/${COHORT}.saf.idx"
    if [[ ! -f "$SAF" ]]; then
        echo "WARNING: $SAF not found, skipping $COHORT (run pop_SAF.sh first)" >&2
        continue
    fi

    echo "[$(date)] Starting $COHORT"

    # 1. Folded global SFS (same realSFS call as pop_het.sh)
    SFS="${OUTDIR}/${COHORT}.folded.sfs"
    if [[ -s "$SFS" ]]; then
        echo "[$(date)] $COHORT: folded SFS exists, skipping realSFS"
    else
        echo "[$(date)] Step 1: realSFS for $COHORT"
        "$RS" "$SAF" -fold 1 -P "$CORES" > "$SFS"
    fi

    # 2. Per-site thetas from the SAF + folded SFS prior
    echo "[$(date)] Step 2: saf2theta for $COHORT"
    "$MISC/realSFS" saf2theta "$SAF" -sfs "$SFS" -fold 1 -outname "${OUTDIR}/${COHORT}" 2> "${OUTDIR}/logs/${COHORT}.saf2theta.log"

    # 3. Genome-wide (per-contig) summary stats, incl. Watterson's theta
    echo "[$(date)] Step 3: thetaStat do_stat for $COHORT"
    "$MISC/thetaStat" do_stat "${OUTDIR}/${COHORT}.thetas.idx" 2> "${OUTDIR}/logs/${COHORT}.thetaStat.log"

    # 4. Summarize genome-wide Watterson's theta -> Ne
    #    pestPG columns (no header token on col1): ... tW tP tF tH tL Tajima fuf fud fayh zeng nSites
    #    tW = per-window sum of Watterson's theta; nSites = sites with data in window.
    PESTPG="${OUTDIR}/${COHORT}.thetas.idx.pestPG"
    read sum_tW sum_nSites <<< "$(awk 'NR>1 {tw+=$4; ns+=$14} END{print tw, ns}' "$PESTPG")"
    theta_w_per_site=$(awk -v tw="$sum_tW" -v ns="$sum_nSites" 'BEGIN{printf "%.10g", tw/ns}')
    ne=$(awk -v tw="$theta_w_per_site" -v mu="$MU" 'BEGIN{printf "%.2f", tw/(4*mu)}')
    echo -e "${COHORT}\t${sum_tW}\t${sum_nSites}\t${theta_w_per_site}\t${ne}" >> "$SUMMARY"

    echo "[$(date)] Done: $COHORT"
done

echo "Wrote summary to $SUMMARY"
cat "$SUMMARY"
