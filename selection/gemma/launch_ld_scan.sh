#!/bin/bash
# Submits ld_scan_array.sh over scaffolds carrying a candidate gene
# (candidate_gene_intervals.tsv from find_intergenic_candidate_genes.py),
# split into small/large jobs by LD file size.
# Usage: bash launch_ld_scan.sh [r2_threshold, default 0.7]  (run with bash, not sbatch)

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

LD_DIR=${BASE_DIR}/variant_calling/ngsLD
WORKDIR=${BASE_DIR}/selection/gemma/ld_curation
ARRAY_SCRIPT=${BASE_DIR}/gvtp_popgen/selection/gemma/ld_scan_array.sh
R2_THRESHOLD="${1:-0.7}"

cd "$WORKDIR"
mkdir -p linked_positions logs

if [[ ! -f candidate_gene_intervals.tsv ]]; then
    echo "candidate_gene_intervals.tsv not found -- run find_intergenic_candidate_genes.py first" >&2
    exit 1
fi

cut -f1 candidate_gene_intervals.tsv | sort -u > scaffolds_needed.txt

> scaffolds_small.txt
> scaffolds_large.txt
while read -r scaffold; do
    ld_gz="${LD_DIR}/${scaffold}_prune.ld.gz"
    if [[ ! -f "$ld_gz" ]]; then
        echo "WARNING: no LD file for ${scaffold}, skipping" >&2
        continue
    fi
    if [[ -f "linked_positions/${scaffold}.linked.tsv" ]]; then
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
