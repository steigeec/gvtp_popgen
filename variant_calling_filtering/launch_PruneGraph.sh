#!/bin/bash
# Submit PruneGraph_array.sh array jobs for all ngsLD outputs, split by file size.
# Usage: bash launch_PruneGraph.sh  (not sbatch; the array size is set at submission)

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

WORKDIR=${BASE_DIR}/variant_calling/ngsLD
cd "$WORKDIR"
mkdir -p pos logs

find . -name "*_prune.ld" -not -size +1G | sed 's|^./||' | sort > ld_small_filelist.txt
find . -name "*_prune.ld"      -size +1G | sed 's|^./||' | sort > ld_large_filelist.txt

N_SMALL=$(wc -l < ld_small_filelist.txt)
N_LARGE=$(wc -l < ld_large_filelist.txt)
echo "Small (<=1 GB): $N_SMALL files"
echo "Large  (>1 GB): $N_LARGE files"

# small files: up to 200 concurrent tasks
if (( N_SMALL > 0 )); then
    JID_SMALL=$(sbatch --parsable \
        --array=1-${N_SMALL}%200 \
        --mem=5G \
        --time=2:00:00 \
        PruneGraph_array.sh ld_small_filelist.txt)
    echo "Submitted small array job: $JID_SMALL"
fi

# large files: up to 4 concurrent tasks
if (( N_LARGE > 0 )); then
    JID_LARGE=$(sbatch --parsable \
        --array=1-${N_LARGE}%4 \
        --mem=30G \
        --time=24:00:00 \
        PruneGraph_array.sh ld_large_filelist.txt)
    echo "Submitted large array job: $JID_LARGE"
fi
