#!/bin/bash
# Split sites not already paralog-filtered into chunks, submit ngsParalog_chunk.sh
# as an array job, then merge the per-chunk .lr files.
# Usage: bash launch_ngsParalog.sh  (submits the jobs itself)

BASE_DIR=/private/groups/meyerlab/esteiger

WORKDIR="${BASE_DIR}/variant_calling/ngsParalog"
SNPS_DIR="${BASE_DIR}/variant_calling/angsd"
MAFS="${SNPS_DIR}/noMAF_preParalog.mafs.gz"
PREV_PARALOG="${BASE_DIR}/variant_calling/angsd/paralog_filtered/postNGSparalog.sites.gz"
CHUNKS_DIR="${WORKDIR}/chunks"
N_CHUNKS=30

mkdir -p "$CHUNKS_DIR" "${WORKDIR}/logs"

# ---------------------------------------------------------------------------
echo "=== Step 1: Extracting sites, separating known-good from new ==="
zcat "$MAFS" | grep -v "^chromo" | cut -f1-2 > "${WORKDIR}/all_sites.txt"
NSITES=$(wc -l < "${WORKDIR}/all_sites.txt")
echo "  Total sites: $NSITES"

# sites already in postNGSparalog.sites.gz passed ngsParalog before and are not re-run
PREV_FILTERED="${WORKDIR}/previously-paralog-filtered.sites"
NEW_SITES="${WORKDIR}/new_sites_tmp.txt"
awk -v prev_file="$PREV_FILTERED" -v new_file="$NEW_SITES" \
    'NR==FNR { seen[$1"\t"$2]=1; next }
     ($1"\t"$2) in seen { print > prev_file; next }
     { print > new_file }' \
  <(zcat "$PREV_PARALOG") "${WORKDIR}/all_sites.txt"
echo "  Previously filtered (skipping): $(wc -l < "$PREV_FILTERED")"
echo "  New sites to evaluate: $(wc -l < "$NEW_SITES")"

split -n l/${N_CHUNKS} "$NEW_SITES" "${CHUNKS_DIR}/sites_chunk_"
rm "$NEW_SITES"
ACTUAL_CHUNKS=$(ls "${CHUNKS_DIR}"/sites_chunk_* | wc -l)
echo "  Chunks created: $ACTUAL_CHUNKS"

# ---------------------------------------------------------------------------
echo ""
echo "=== Step 2: Submitting array job (${ACTUAL_CHUNKS} tasks) ==="
ARRAY_JOB=$(sbatch --parsable \
  --array=0-$((ACTUAL_CHUNKS - 1)) \
  "${WORKDIR}/ngsParalog_chunk.sh" "${CHUNKS_DIR}")
echo "  Array job ID: ${ARRAY_JOB}"

# ---------------------------------------------------------------------------
echo ""
echo "=== Step 3: Submitting merge job (runs after all chunks succeed) ==="
MERGE_JOB=$(sbatch --parsable \
  --dependency=afterok:${ARRAY_JOB} \
  --job-name=ngsParalog_merge \
  --partition=long \
  --ntasks=1 --mem=16gb --time=2:00:00 \
  --output="${WORKDIR}/logs/merge_%j.log" \
  --wrap="cat '${CHUNKS_DIR}'/sites_chunk_*.lr \
    | sort -k1,1 -k2,2n > '${WORKDIR}/noMAF_preParalog.lr' \
    && echo \"Total sites: \$(wc -l < '${WORKDIR}/noMAF_preParalog.lr')\"")
echo "  Merge job ID: ${MERGE_JOB}"

echo ""
echo "=== Submitted. Monitor with: squeue -u \$USER ==="
