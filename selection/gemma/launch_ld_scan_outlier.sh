#!/bin/bash
# Run as: bash launch_ld_scan_outlier.sh  (NOT sbatch -- array size must be known at submission)
#
# Submits ld_scan_outlier_array.sh over only the scaffolds that actually
# carry a fog/BIO18/BIO15 outlier SNP (outlier_snp_seeds.tsv) -- not the
# whole genome. Split into small/large tiers by LD file size, same pattern
# as launch_ld_scan.sh (the candidate-gene LD scan), since it's already
# proven at this scale on these exact files.

set -euo pipefail

LD_DIR=/private/groups/meyerlab/esteiger/variant_calling/ngsLD
WORKDIR=/private/groups/meyerlab/esteiger/selection/gemma/ld_outlier_windows
ARRAY_SCRIPT=/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/ld_scan_outlier_array.sh
R2_THRESHOLD="${1:-0.4}"

cd "$WORKDIR"
mkdir -p linked_snp_positions logs

if [[ ! -f outlier_snp_seeds.tsv ]]; then
    echo "outlier_snp_seeds.tsv not found -- build it first (var_name, scaffold, pos per outlier SNP)" >&2
    exit 1
fi

cut -f2 outlier_snp_seeds.tsv | sort -u > scaffolds_needed.txt

> scaffolds_small.txt
> scaffolds_large.txt
while read -r scaffold; do
    ld_gz="${LD_DIR}/${scaffold}_prune.ld.gz"
    if [[ ! -f "$ld_gz" ]]; then
        echo "WARNING: no LD file for ${scaffold}, skipping" >&2
        continue
    fi
    if [[ -f "linked_snp_positions/${scaffold}.linked.tsv" ]]; then
        continue
    fi
    size=$(stat -c%s "$ld_gz")
    if (( size > 1073741824 )); then
        echo "$scaffold" >> scaffolds_large.txt
    else
        echo "$scaffold" >> scaffolds_small.txt
    fi
done < scaffolds_needed.txt

N_SMALL=$(wc -l < scaffolds_small.txt)
N_LARGE=$(wc -l < scaffolds_large.txt)
echo "Small (<=1 GB): $N_SMALL scaffolds"
echo "Large  (>1 GB): $N_LARGE scaffolds"

if (( N_SMALL > 0 )); then
    JID_SMALL=$(sbatch --parsable \
        --array=1-${N_SMALL}%20 \
        --cpus-per-task=4 \
        --mem=8G \
        --time=2:00:00 \
        "$ARRAY_SCRIPT" scaffolds_small.txt "$R2_THRESHOLD")
    echo "Submitted small array job: $JID_SMALL"
fi

if (( N_LARGE > 0 )); then
    JID_LARGE=$(sbatch --parsable \
        --array=1-${N_LARGE}%4 \
        --cpus-per-task=8 \
        --mem=16G \
        --time=24:00:00 \
        "$ARRAY_SCRIPT" scaffolds_large.txt "$R2_THRESHOLD")
    echo "Submitted large array job: $JID_LARGE"
fi
