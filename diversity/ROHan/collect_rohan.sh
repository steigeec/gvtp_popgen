#!/bin/bash
# collect_rohan.sh
# Concatenate per-sample ROHan output into summary tables. Run after all
# rohan_phase1.sh/rohan_phase2.sh tasks have finished.
#
# ROHan pass 1 (MCMC) gives theta as a point estimate with lower/upper bounds
# (thetamid thetamin thetamax in each .summary.txt), and pass 2 (HMM) was run
# once per value:
#   .mid.hmmrohl.gz   theta point estimate
#   .min.hmmrohl.gz   lower bound of theta, so more windows fall below
#                     background and more ROH are called
#   .max.hmmrohl.gz   upper bound of theta, so fewer ROH are called
# The spread reflects uncertainty in background heterozygosity (mostly driven
# by data volume per sample), not confidence in the ROH calls, and the ROHan
# README notes these bounds tend to be inflated. .mid is used as the primary
# call set; .min/.max are kept only as per-sample totals for spread plots.
# Samples with no ROH contain the literal string "none found".
#
# Outputs (in OUTDIR):
#   rohan_regions.txt                   mid-theta segments, input for roh_plot.R
#   rohan_summary_stats.txt             theta estimates and % genome in ROH
#   rohan_roh_burden_by_threshold.txt   total ROH Mb per sample at mid/min/max

#SBATCH --job-name=ROHan_collect
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=1:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/rohan_collect_%j.log

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

OUTDIR=${BASE_DIR}/diversity/ROH/ROHan

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

# Skip the header line (#ROH_ID CHROM BEGIN END ROH_LENGTH VALIDATED_SITES)
    echo "$content" \
        | awk -v s="$sample" 'NR>1 {OFS="\t"; print s, $2, $3, $4, $5, $6}' \
        >> "$OUTDIR/rohan_regions.txt"
    (( n_with_roh++ )) || true
done

n_segs_raw=$(( $(wc -l < "$OUTDIR/rohan_regions.txt") - 1 ))
echo "ROH segments (raw): $n_segs_raw across $n_with_roh samples ($n_no_roh samples had no ROH)"

# --- Merge segments separated by at most one 100 kb window --------------------
# ROHan already merges contiguous windows, but one noisy window dropping out of
# the ROH state can split a run into two segments separated by exactly one
# --size window. Merging across a 1-window gap fixes that (~2% fewer segments);
# anything wider is left alone.
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
