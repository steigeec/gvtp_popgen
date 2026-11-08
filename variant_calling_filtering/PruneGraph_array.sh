#!/bin/bash

#SBATCH --job-name=PruneGraph
#SBATCH --partition=long
#SBATCH --cpus-per-task=40
#SBATCH --output=logs/PruneGraph_%A_%a.log

# prune_graph on one ngsLD output (array task): drop linked sites
# (r2 >= 0.4 within 50 kb). Scaffolds with no linked pairs keep all sites.
# Usage: submitted by launch_PruneGraph.sh

BASE_DIR=/private/groups/meyerlab/esteiger

export TMPDIR=${BASE_DIR}/tmp/

WORKDIR=${BASE_DIR}/variant_calling/ngsLD
cd "$WORKDIR"
mkdir -p pos

PG=${BASE_DIR}/software/prune_graph/target/release/prune_graph
SITES_DB=${BASE_DIR}/variant_calling/angsd/postParalog_MAF.sites

FILELIST=$1
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
NAME="${FILE%_prune.ld}"

if [[ -z "$FILE" ]]; then
    echo "No file for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

if [[ -f "pos/${NAME}_prune.pos" ]]; then
    echo "Already done: ${NAME}"
    exit 0
fi

if [[ ! -f "$FILE" ]]; then
    echo "File not found (already gzipped?): $FILE"
    exit 0
fi

pg_stderr=$($PG --in "$FILE" -n 40 --header \
    --weight-field "r2" \
    --weight-filter "dist <= 50000 && r2 >= 0.4" \
    --out "pos/${NAME}_prune.pos" 2>&1)
pg_exit=$?

if [[ $pg_exit -eq 0 ]]; then
    pigz -p 40 "$FILE"
    echo "Done: ${NAME}"
elif echo "$pg_stderr" | grep -q "Graph is empty"; then
    grep "^${NAME}	" "$SITES_DB" \
        | awk '{print $1":"$2}' > "pos/${NAME}_prune.pos"
    pigz -p 40 "$FILE"
    echo "Graph is empty for ${NAME} — wrote all sites"
else
    echo "ERROR for ${NAME}: $pg_stderr"
    exit 1
fi
