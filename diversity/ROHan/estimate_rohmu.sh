#!/bin/bash
# estimate_rohmu.sh
# Calibrate ROHan's rohmu parameter. Selects the 2 highest-coverage samples per
# cluster (1 for Dangermond) so the calibration set spans the expected range of
# inbreeding, builds 3x subsampled BAMs from the 6x BAMs, writes the phase 1
# task list, and submits rohmu_phase1.sh and rohmu_phase2.sh.
#
# Requires: prep_samples.sh (subsampled_bams/)
# Usage:    bash estimate_rohmu.sh
#           srun --mem=4gb --time=1:00:00 bash estimate_rohmu.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

FIELD_DATA=${BASE_DIR}/field_data.txt
OUTDIR=${BASE_DIR}/diversity/ROH/ROHan
SCRIPT_DIR=${BASE_DIR}/gvtp_popgen/diversity/ROHan
SUBSAMPLE_DIR=$OUTDIR/subsampled_bams
ROHMU_DIR=$OUTDIR/rohmu_calibration

CLUSTERS=(CDFW VAFB_Lion_Head VAFB_Arguello Dangermond Point_Conception Strauss)
SAMPLES_PER_CLUSTER=2
MIN_COV=6
CAL_CHAINS=5000
CAL_SEED=42

mkdir -p "$ROHMU_DIR/logs"

# --- Stratified sample selection -----------------------------------------------
> "$ROHMU_DIR/cal_samples.tsv"
for cluster in "${CLUSTERS[@]}"; do
    n=$SAMPLES_PER_CLUSTER
    [[ "$cluster" == "Dangermond" ]] && n=1
    awk -F'\t' -v cl="$cluster" -v min_cov="$MIN_COV" \
        -v sub_dir="$SUBSAMPLE_DIR" '
    NR==1 {
        for (i=1; i<=NF; i++) {
            if ($i == "SampleID") sid_col = i
            if ($i == "Cov")      cov_col = i
            if ($i == "Clusters") clust_col = i
        }
        next
    }
    $sid_col == "4531" || $sid_col == "7930" || $sid_col == "7976" { next }
    $clust_col == cl && $cov_col+0 >= min_cov {
        printf "%.4f\t%s\t%s/%s_sub6x.bam\n",
            $cov_col+0, $sid_col, sub_dir, $sid_col
    }
    ' "$FIELD_DATA" | sort -rn | head -"$n"
done >> "$ROHMU_DIR/cal_samples.tsv"

n_cal=$(wc -l < "$ROHMU_DIR/cal_samples.tsv")
echo "Calibration samples: $n_cal selected (stratified by cluster)"
cat "$ROHMU_DIR/cal_samples.tsv"
echo ""

if (( n_cal == 0 )); then
    echo "ERROR: no calibration samples selected" >&2; exit 1
fi

# --- Verify 6x BAMs exist -----------------------------------------------------
while IFS=$'\t' read -r cov sample sub6x_bam; do
    [[ -s "$sub6x_bam" ]] || {
        echo "ERROR: 6x BAM missing for $sample ($sub6x_bam). Run prep_samples.sh first." >&2
        exit 1
    }
done < "$ROHMU_DIR/cal_samples.tsv"

# --- Build 3x subsampled BAMs -------------------------------------------------
while IFS=$'\t' read -r cov sample sub6x_bam; do
    local_sub3x="$ROHMU_DIR/${sample}_sub3x.bam"
    if [[ ! -s "$local_sub3x" ]]; then
        echo "  Building 3x BAM for $sample ..."
        samtools view -@ 1 -s "${CAL_SEED}.500000" -b -o "$local_sub3x" "$sub6x_bam"
        samtools index "$local_sub3x"
    fi
done < "$ROHMU_DIR/cal_samples.tsv"
echo "3x BAMs ready."

# --- Write phase 1 task list --------------------------------------------------
{
    while IFS=$'\t' read -r cov sample sub6x_bam; do
        printf '%s\t%s\t%d\n' \
            "${ROHMU_DIR}/${sample}_6x" "$sub6x_bam" "$CAL_CHAINS"
        printf '%s\t%s\t%d\n' \
            "${ROHMU_DIR}/${sample}_3x" "${ROHMU_DIR}/${sample}_sub3x.bam" "$CAL_CHAINS"
    done < "$ROHMU_DIR/cal_samples.tsv"
} > "$ROHMU_DIR/phase1_tasks.tsv"

n_tasks=$(wc -l < "$ROHMU_DIR/phase1_tasks.tsv")
echo "Phase 1 tasks: $n_tasks (written to phase1_tasks.tsv)"

# --- Submit SLURM jobs ---------------------------------------------------------
echo ""
echo "Submitting phase 1 array ($n_tasks tasks)..."
P1=$(sbatch --parsable --array=1-"${n_tasks}" "$SCRIPT_DIR/rohmu_phase1.sh")
echo "  Phase 1 job ID: $P1"

echo "Submitting phase 2 (depends on phase 1)..."
P2=$(sbatch --parsable --dependency=afterok:"$P1" "$SCRIPT_DIR/rohmu_phase2.sh")
echo "  Phase 2 job ID: $P2"

echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Phase 2 runs automatically after all phase 1 tasks complete."
