#!/bin/bash
# prep_samples.sh
# Select samples with coverage > 3x (excluding outlier 4531) and downsample BAMs
# above 6x to 6x (fixed seed). BAMs between 3x and 6x are used as is.
#
# Outputs (in OUTDIR):
#   samples_gt3x.tsv    sample, original coverage, original BAM, analysis BAM
#   bamlist_rohan.txt   analysis BAM paths
#   subsampled_bams/    downsampled BAMs and indices
#
# Run before estimate_rohmu.sh and launch_rohan.sh.

#SBATCH --job-name=prep_rohan
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16gb
#SBATCH --time=168:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/prep_rohan_%j.log

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

BAMDIR=${BASE_DIR}/mapping/merged
FIELD_DATA=${BASE_DIR}/field_data.txt
OUTDIR=${BASE_DIR}/diversity/ROH/ROHan
FAI=${BASE_DIR}/genome/DEINVI_genome.fasta.masked.fai
SUBSAMPLE_DIR=$OUTDIR/subsampled_bams
TARGET_COV=6
SEED=42
N_PARALLEL=${SLURM_NTASKS:-4}

mkdir -p "$SUBSAMPLE_DIR"

# Scaffolds >= 1 Mb
awk '$2 >= 1000000 {print $1}' "$FAI" > "$OUTDIR/scaffolds_1mb_names.txt"

# --- Build sample table: cov > 3x, excluding 4531 ----------------------------
awk -v bamdir="$BAMDIR" -v sub_dir="$SUBSAMPLE_DIR" -v target="$TARGET_COV" '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "SampleID") sid_col = i
        if ($i == "Cov")      cov_col = i
    }
    next
}
$sid_col == "4531" { next }
$cov_col + 0 > 3 {
    sample = $sid_col
    cov    = $cov_col + 0
    orig   = bamdir "/" sample "_rmdup.bam"
    final  = (cov > target) ? sub_dir "/" sample "_sub" target "x.bam" : orig
    printf "%s\t%.4f\t%s\t%s\n", sample, cov, orig, final
}
' "$FIELD_DATA" | sort -t$'\t' -k2 -rn > "$OUTDIR/samples_gt3x.tsv"

n_total=$(wc -l < "$OUTDIR/samples_gt3x.tsv")
n_sub=$(awk -v t="$TARGET_COV" '$2 > t' "$OUTDIR/samples_gt3x.tsv" | wc -l)
echo "Samples passing >3x filter: $n_total ($n_sub will be subsampled to ${TARGET_COV}x)"

# --- Per-sample: subsample or index as needed --------------------------------
process_sample() {
    local sample=$1 cov=$2 orig_bam=$3 final_bam=$4

    if [[ ! -f "$orig_bam" ]]; then
        echo "WARNING: BAM not found for $sample — skipping" >&2
        return 0
    fi

    if [[ "$final_bam" != "$orig_bam" ]]; then
        local frac s_arg
        frac=$(awk -v c="$cov" -v t="$TARGET_COV" 'BEGIN {printf "%.6f", t/c}')
        s_arg="${SEED}${frac#0}"   # samtools format: integer_seed.fraction (e.g. 42.183275)
        if [[ ! -s "$final_bam" ]]; then
            echo "  Subsampling $sample (${cov}x → ${TARGET_COV}x, -s ${s_arg}) ..." >&2
            samtools view -@ 1 -s "$s_arg" -b -o "$final_bam" "$orig_bam"
        fi
    fi

    [[ -s "${final_bam}.bai" ]] || samtools index "$final_bam"
    echo "$final_bam"
}

export -f process_sample
export SUBSAMPLE_DIR TARGET_COV SEED

parallel --keep-order -j "$N_PARALLEL" --colsep '\t' \
    'process_sample {1} {2} {3} {4}' \
    :::: "$OUTDIR/samples_gt3x.tsv" \
    > "$OUTDIR/bamlist_rohan.txt"

echo "bamlist_rohan.txt written: $(wc -l < "$OUTDIR/bamlist_rohan.txt") analysis-ready BAMs"
echo "Next: sbatch estimate_rohmu.sh  (then review chosen_rohmu.txt before launching full run)"
