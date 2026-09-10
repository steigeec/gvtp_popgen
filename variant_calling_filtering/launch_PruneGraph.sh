#!/bin/bash
# Run as: bash launch_PruneGraph.sh  (NOT sbatch -- array size must be known at submission)

set -euo pipefail

WORKDIR=/private/groups/meyerlab/esteiger/variant_calling/ngsLD
cd "$WORKDIR"
mkdir -p pos logs

# Split unprocessed .ld files into two tiers by size
find . -name "*_prune.ld" -not -size +1G | sed 's|^./||' | sort > ld_small_filelist.txt
find . -name "*_prune.ld"      -size +1G | sed 's|^./||' | sort > ld_large_filelist.txt

N_SMALL=$(wc -l < ld_small_filelist.txt)
N_LARGE=$(wc -l < ld_large_filelist.txt)
echo "Small (<=1 GB): $N_SMALL files"
echo "Large  (>1 GB): $N_LARGE files"

# Small files: up to 200 concurrent tasks, 30 GB each, short time limit
if (( N_SMALL > 0 )); then
    JID_SMALL=$(sbatch --parsable \
        --array=1-${N_SMALL}%200 \
        --mem=5G \
        --time=2:00:00 \
        PruneGraph_array.sh ld_small_filelist.txt)
    echo "Submitted small array job: $JID_SMALL"
fi

# Large files: up to 4 concurrent tasks, 300 GB each, long time limit
if (( N_LARGE > 0 )); then
    JID_LARGE=$(sbatch --parsable \
        --array=1-${N_LARGE}%4 \
        --mem=30G \
        --time=24:00:00 \
        PruneGraph_array.sh ld_large_filelist.txt)
    echo "Submitted large array job: $JID_LARGE"
fi
