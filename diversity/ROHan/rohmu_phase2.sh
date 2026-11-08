#!/bin/bash
# rohmu_phase2.sh
# ROHan pass 2 (HMM) for each calibration sample at 6x and 3x across candidate
# rohmu values, then pick the rohmu with the lowest mean relative cross-depth
# discordance in F_ROH (|dF_ROH| / mean F_ROH). Writes chosen_rohmu.txt.
#
# Relative rather than absolute discordance is used because F_ROH shrinks at low
# rohmu, so absolute differences are trivially smaller there and always favour
# the lowest value tested. The range extends down to 1e-5 to confirm the minimum
# lies inside the tested range.
#
# Submit via estimate_rohmu.sh or:
#   sbatch --dependency=afterok:<phase1_jobid> rohmu_phase2.sh

#SBATCH --job-name=rohmu_p2
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=50gb
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/rohmu_calibration/logs/phase2_%j.log

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

ROHAN=${BASE_DIR}/software/rohan/bin/rohan
REF=${BASE_DIR}/genome/DEINVI_genome.fasta.masked
OUTDIR=${BASE_DIR}/diversity/ROH/ROHan
ROHMU_DIR=$OUTDIR/rohmu_calibration
SCAFFOLDS=$OUTDIR/scaffolds_1mb_names.txt

ROHMU_VALUES=(1e-5 2e-5 3e-5 5e-5 8e-5 2e-4 5e-4 8e-4 2e-3)
CAL_CHAINS=5000
N_PARALLEL=${SLURM_NTASKS:-10}

# --- Pass 2 helper ------------------------------------------------------------
run_pass2() {
    local base_prefix=$1 rohmu=$2 bam=$3
    local out_prefix="${base_prefix}_rohmu${rohmu}"
    local hest_target
    hest_target=$(basename "${base_prefix}.hEst.gz")
    [[ -e "${out_prefix}.hEst.gz" ]] || \
        ln -s "$hest_target" "${out_prefix}.hEst.gz"
    if [[ ! -s "${out_prefix}.summary.txt" ]]; then
        "$ROHAN" \
            --tstv 1.86 --size 100000 --chains "$CAL_CHAINS" \
            --hmm --rohmu "$rohmu" \
            --auto "$SCAFFOLDS" \
            -o "$out_prefix" "$REF" "$bam" \
            2>> "${out_prefix}.log"
    fi
}
export -f run_pass2
export ROHAN REF SCAFFOLDS CAL_CHAINS

# --- Phase 2: HMM ROH calling across all rohmu values -------------------------
echo "=== Phase 2: HMM calling across rohmu values ==="

declare -a pass2_jobs=()
while IFS=$'\t' read -r cov sample sub6x_bam; do
    for rohmu in "${ROHMU_VALUES[@]}"; do
        pass2_jobs+=("${ROHMU_DIR}/${sample}_6x|${rohmu}|${sub6x_bam}")
        pass2_jobs+=("${ROHMU_DIR}/${sample}_3x|${rohmu}|${ROHMU_DIR}/${sample}_sub3x.bam")
    done
done < "$ROHMU_DIR/cal_samples.tsv"

printf '%s\n' "${pass2_jobs[@]}" \
| parallel -j "$N_PARALLEL" --colsep '\|' \
    --joblog "$ROHMU_DIR/phase2.joblog" \
    'run_pass2 {1} {2} {3}'

echo "Phase 2 complete."

# --- Extract F_ROH from all summary files --------------------------------------
echo "=== Extracting F_ROH ==="

{
    printf 'sample\tdepth\trohmu\tF_ROH\n'
    while IFS=$'\t' read -r cov sample sub6x_bam; do
        for depth in 6x 3x; do
            for rohmu in "${ROHMU_VALUES[@]}"; do
                sumf="${ROHMU_DIR}/${sample}_${depth}_rohmu${rohmu}.summary.txt"
                [[ -s "$sumf" ]] || continue
                froh=$(awk '
                    /Segments in ROH\(%\).*:/ {
                        match($0, /[0-9]+\.[0-9]+/)
                        print substr($0, RSTART, RLENGTH)+0
                    }' "$sumf")
                printf '%s\t%s\t%s\t%s\n' "$sample" "$depth" "$rohmu" "${froh:-NA}"
            done
        done
    done < "$ROHMU_DIR/cal_samples.tsv"
} > "$ROHMU_DIR/froh_by_rohmu.tsv"

echo "F_ROH table → rohmu_calibration/froh_by_rohmu.tsv"

# --- Compute cross-depth discordance and select best rohmu --------------------
# Pairs where both depths give 0% ROH are skipped: they are uninformative, and
# including them would inflate apparent precision at low rohmu, where the HMM
# calls nothing as ROH.
awk -v chosen="$OUTDIR/chosen_rohmu.txt" '
BEGIN { OFS="\t" }
NR==1 { next }
$4 == "NA" { next }
{
    key = $1 SUBSEP $3
    if ($2 == "6x") f6[key] = $4+0
    if ($2 == "3x") f3[key] = $4+0
}
END {
    print "rohmu\tmean_rel_discordance\tn_samples" > "/dev/stderr"
    best_rohmu = ""; best_disc = 1e9
    for (k in f6) {
        if (!(k in f3)) continue
        split(k, a, SUBSEP)
        rohmu = a[2]
        avg = (f6[k] + f3[k]) / 2
        if (avg == 0) continue
        diff = f6[k] - f3[k]; if (diff < 0) diff = -diff
        sum_diff[rohmu] += diff / avg
        n[rohmu]++
    }
    for (rohmu in sum_diff) {
        disc = sum_diff[rohmu] / n[rohmu]
        printf "%s\t%.4f\t%d\n", rohmu, disc, n[rohmu] > "/dev/stderr"
        if (disc < best_disc) { best_disc = disc; best_rohmu = rohmu }
    }
    print best_rohmu > chosen
    printf "\nBest rohmu: %s  (mean |ΔF_ROH|/F_ROH = %.4f)\n", best_rohmu, best_disc > "/dev/stderr"
}
' "$ROHMU_DIR/froh_by_rohmu.tsv" 2>&1 | tee "$ROHMU_DIR/rohmu_discordance.tsv"

echo ""
echo "chosen_rohmu.txt → $(cat "$OUTDIR/chosen_rohmu.txt")"
echo ""
echo "REVIEW BEFORE PROCEEDING:"
echo "  rohmu_calibration/froh_by_rohmu.tsv    — F_ROH per sample/depth/rohmu"
echo "  rohmu_calibration/rohmu_discordance.tsv — relative discordance summary"
echo "Automated selection minimises relative cross-depth discordance (|ΔF_ROH|/F_ROH)."
echo "Confirm F_ROH values are biologically plausible before running launch_rohan.sh."
