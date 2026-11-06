#!/bin/bash
# =============================================================================
# collect_rohan.sh — Concatenate ROHan ROH segment files into rohan_regions.txt
#
# ROHan's pass 1 (MCMC) doesn't produce a single heterozygosity (theta)
# estimate — it produces a credible interval (point estimate + lower/upper
# bounds; see the "thetamid thetamin thetamax" triplet in each .summary.txt).
# Pass 2 (HMM) is then run three separate times, once per bound, each writing
# its own ROH call file:
#   .mid.hmmrohl.gz  — HMM run with the POINT ESTIMATE of theta (best guess)
#   .min.hmmrohl.gz  — HMM run with the LOWER bound of theta -> more windows
#                      look "significantly below background" -> MORE ROH calls
#   .max.hmmrohl.gz  — HMM run with the UPPER bound of theta -> fewer windows
#                      look unusual by comparison -> FEWER ROH calls
# These are not posterior-probability confidence thresholds on the ROH calls
# themselves; the min/max spread reflects uncertainty in the background
# heterozygosity estimate (itself driven mostly by per-sample data volume),
# and ROHan's own README/FAQ notes these bounds tend to run inflated (they're
# computed from the second derivative of a likelihood that's often flat near
# 0). So .mid is used as the primary/canonical call set below; .min/.max are
# additionally collected per-sample (total Mb only, not full segment lists)
# so downstream plots can show the mid estimate with a min/max spread.
#
# Per-segment columns: ROH_ID, CHROM, BEGIN, END, ROH_LENGTH, VALIDATED_SITES
# Samples with no ROH contain the literal string "none found".
#
# Run this once launch_rohan.sh reports no samples remaining (all
# rohan_phase1.sh/rohan_phase2.sh array tasks have completed).
# Output: rohan_regions.txt                    (mid-threshold segments, input for roh_plot.R)
#         rohan_summary_stats.txt              (theta estimates, % genome in ROH per individual)
#         rohan_roh_burden_by_threshold.txt    (per-sample total ROH Mb at mid/min/max, for spread plots)
# =============================================================================

#SBATCH --job-name=ROHan_collect
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=1:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/rohan_collect_%j.log

set -euo pipefail

OUTDIR=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan

# --- Collect ROH segments (mid / point-estimate threshold) -------------------
echo "Collecting ROH segments (mid / point-estimate threshold) ..."
echo -e "Sample\tChrom\tBegin\tEnd\tROH_length_bp\tValidated_sites" \
    > "$OUTDIR/rohan_regions.txt"

n_with_roh=0
n_no_roh=0

for rohl in "$OUTDIR"/*.mid.hmmrohl.gz; do
    [[ -f "$rohl" ]] || continue
    sample=$(basename "$rohl" .mid.hmmrohl.gz)
    [[ "$sample" == test_* ]] && continue

    content=$(zcat "$rohl")
    if [[ "$content" == "none found" ]]; then
        (( n_no_roh++ )) || true
        continue
    fi

    # File has header line (#ROH_ID CHROM BEGIN END ROH_LENGTH VALIDATED_SITES)
    # Skip header, prepend sample name
    echo "$content" \
        | awk -v s="$sample" 'NR>1 {OFS="\t"; print s, $2, $3, $4, $5, $6}' \
        >> "$OUTDIR/rohan_regions.txt"
    (( n_with_roh++ )) || true
done

n_segs_raw=$(( $(wc -l < "$OUTDIR/rohan_regions.txt") - 1 ))
echo "ROH segments (raw): $n_segs_raw across $n_with_roh samples ($n_no_roh samples had no ROH)"

# --- Merge near-adjacent segments (same sample+scaffold, gap <= 1 window) ----
# ROHan's own HMM output already merges truly-contiguous windows, but a single
# spurious window flipping out of the ROH state (e.g. one noisy heterozygous
# call) can fragment what's otherwise one real run into two reported segments
# separated by exactly one 100kb (--size) window. Merging across a 1-window
# gap is a conservative fix for that specific artifact — it barely changes
# the overall length distribution (checked empirically: ~2% fewer segments)
# but avoids undercounting genuinely-contiguous long runs. NOT intended to
# paper over genuine short separate ROH — anything wider than 1 window is
# left alone.
MERGE_GAP_TOL_BP=100000
python3 - "$OUTDIR/rohan_regions.txt" "$MERGE_GAP_TOL_BP" << 'PYEOF'
import sys
from collections import defaultdict

path, tol = sys.argv[1], int(sys.argv[2])
with open(path) as f:
    header = f.readline()
    rows = [line.rstrip("\n").split("\t") for line in f]

groups = defaultdict(list)
for sample, chrom, begin, end, length, valid in rows:
    groups[(sample, chrom)].append((int(begin), int(end), int(valid)))

merged = []
for (sample, chrom), ivs in groups.items():
    ivs.sort()
    cur_b, cur_e, cur_v = ivs[0]
    for b, e, v in ivs[1:]:
        if b - cur_e - 1 <= tol:
            cur_e = e
            cur_v += v
        else:
            merged.append((sample, chrom, cur_b, cur_e, cur_v))
            cur_b, cur_e, cur_v = b, e, v
    merged.append((sample, chrom, cur_b, cur_e, cur_v))

merged.sort(key=lambda r: (r[0], r[1], r[2]))
with open(path, "w") as f:
    f.write(header)
    for sample, chrom, b, e, v in merged:
        f.write(f"{sample}\t{chrom}\t{b}\t{e}\t{e - b + 1}\t{v}\n")
PYEOF

n_segs=$(( $(wc -l < "$OUTDIR/rohan_regions.txt") - 1 ))
echo "ROH segments (after 1-window-gap merge): $n_segs (was $n_segs_raw)"

# --- Collect per-sample total ROH burden at all three thresholds -------------
echo "Collecting total ROH burden at mid/min/max thresholds ..."
echo -e "Sample\tTotal_ROH_Mb_mid\tTotal_ROH_Mb_min\tTotal_ROH_Mb_max" \
    > "$OUTDIR/rohan_roh_burden_by_threshold.txt"

n_burden=0
for midf in "$OUTDIR"/*.mid.hmmrohl.gz; do
    [[ -f "$midf" ]] || continue
    sample=$(basename "$midf" .mid.hmmrohl.gz)
    [[ "$sample" == test_* ]] && continue

    mb_mid=0; mb_min=0; mb_max=0
    for tier in mid min max; do
        f="$OUTDIR/${sample}.${tier}.hmmrohl.gz"
        [[ -f "$f" ]] || continue
        c=$(zcat "$f")
        [[ "$c" == "none found" ]] && continue
        mb=$(echo "$c" | awk 'NR>1 {sum+=$5} END {printf "%.6f", sum/1e6}')
        case "$tier" in
            mid) mb_mid=$mb ;;
            min) mb_min=$mb ;;
            max) mb_max=$mb ;;
        esac
    done
    printf '%s\t%s\t%s\t%s\n' "$sample" "$mb_mid" "$mb_min" "$mb_max" \
        >> "$OUTDIR/rohan_roh_burden_by_threshold.txt"
    (( n_burden++ )) || true
done
echo "Threshold burden collected for $n_burden samples -> rohan_roh_burden_by_threshold.txt"

# --- Collect genome-wide summary stats ---------------------------------------
echo "Collecting genome-wide summary stats ..."
echo -e "Sample\tTheta_outside_ROH\tTheta_outside_lower\tTheta_outside_upper\tTheta_inc_ROH\tTheta_inc_lower\tTheta_inc_upper\tPct_genome_ROH" \
    > "$OUTDIR/rohan_summary_stats.txt"

for sumf in "$OUTDIR"/*.summary.txt; do
    [[ -f "$sumf" ]] || continue
    sample=$(basename "$sumf" .summary.txt)
    [[ "$sample" == test_* ]] && continue

    awk -v s="$sample" '
        /Genome-wide theta outside ROH:/ { split($0,a,"\t"); theta_out=a[2]; lo=a[3]; hi=a[4] }
        /Genome-wide theta inc. ROH:/    { split($0,a,"\t"); theta_inc=a[2]; ilo=a[3]; ihi=a[4] }
        /Segments in ROH\(%\)/           { split($0,a,":"); pct=a[2]; sub(/^[ \t]+/,"",pct); split(pct,b," "); roh_pct=b[1] }
        END { OFS="\t"; print s, theta_out, lo, hi, theta_inc, ilo, ihi, roh_pct }
    ' "$sumf" >> "$OUTDIR/rohan_summary_stats.txt"
done

n_stats=$(( $(wc -l < "$OUTDIR/rohan_summary_stats.txt") - 1 ))
echo "Summary stats: $n_stats samples → rohan_summary_stats.txt"
echo "Done."
