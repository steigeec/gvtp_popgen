#!/bin/bash

#SBATCH --job-name=PruneGraph
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=esteiger@ucsc.edu
#SBATCH --cpus-per-task=40
#SBATCH --output=logs/PruneGraph_%A_%a.log

export TMPDIR=/private/groups/meyerlab/esteiger/tmp/

WORKDIR=/private/groups/meyerlab/esteiger/variant_calling/ngsLD
cd "$WORKDIR"
mkdir -p pos

PG=/private/groups/meyerlab/esteiger/software/prune_graph/target/release/prune_graph
SITES_DB=/private/groups/meyerlab/esteiger/variant_calling/angsd/postParalog_MAF.sites

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
