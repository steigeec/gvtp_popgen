#!/bin/bash
set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

# launch_ngsF.sh
# Run after doGlf_launch.sh. Submits ngsF.sh (100 EM replicates) for each chunk
# not yet complete, then a finalize step (afterany) that picks each chunk's best
# replicate by final logLkl, writes best_LL.txt and ngsF_results.txt, and runs
# inbreeding_plot.R. If nothing needs submitting, finalize runs immediately.
# Usage: bash launch_ngsF.sh [--finalize-only]

NGSF_DIR="${BASE_DIR}/diversity/ngsF"
SCRIPT_DIR="${BASE_DIR}/gvtp_popgen/diversity/inbreeding"

cd "$NGSF_DIR"

finalize () {
  echo "=== Finalize: computing best_LL.txt and ngsF_results.txt ==="
  : > best_LL.txt
  printf 'BAM\tNGSF\n' > ngsF_results.txt

  for DIR in bamlist_chunk_*_output; do
    for file in "$DIR"/ngsF_log_*; do
      final_line=$(grep "Final logLkl: " "$file" | tail -1)
      if [ -n "$final_line" ]; then
        logLkl=$(echo "$final_line" | awk '{print $NF}')
        echo -e "$file\t$logLkl"
      fi
    done > "$DIR/all_ngsF_final_logLkl.txt"

    {
      echo "For $DIR, the best output was:"
      BEST_LINE=$(sort -k2 -nr "$DIR/all_ngsF_final_logLkl.txt" | head -1)
      echo "$BEST_LINE"
    } >> best_LL.txt

    BEST_LOG=$(echo "$BEST_LINE" | cut -f1)
    NUM=$(basename "$BEST_LOG" | sed 's/^ngsF_log_//')
    BAMLIST="${DIR%_output}"
    OUTFILE="${DIR}/500k_sites_${NUM}"
    paste "$BAMLIST" "$OUTFILE" >> ngsF_results.txt
  done

  echo "best_LL.txt and ngsF_results.txt written."

  echo "=== Finalize: running inbreeding_plot.R ==="
  Rscript "${SCRIPT_DIR}/inbreeding_plot.R"
}

if [ "${1:-}" == "--finalize-only" ]; then
  finalize
  exit 0
fi

# --- Stage 3: launch ngsF replicate jobs per chunk (skip chunks already complete) ---
declare -A ngsf_jobid

for chunk in bamlist_chunk_*[!.]*; do
  [[ "$chunk" == *.* ]] && continue

  glf="${chunk}.glf.gz"
  if [ ! -f "$glf" ]; then
    echo "ERROR: $glf not found. Run doGlf_launch.sh first." >&2
    exit 1
  fi

  outdir="${chunk}_output"
  n_done=0
  if [ -d "$outdir" ]; then
    n_done=$(ls "${outdir}"/500k_sites_[0-9]*.pars 2>/dev/null | wc -l)
  fi

  if [ "$n_done" -eq 100 ]; then
    echo "Skipping ngsF for $chunk: 100/100 replicates already complete"
    continue
  fi

  jid=$(sbatch --parsable "${SCRIPT_DIR}/ngsF.sh" "$chunk")
  ngsf_jobid[$chunk]="$jid"
  echo "Submitted ngsF for $chunk (job $jid)"
done

# --- Stage 4: finalize (best_LL.txt, ngsF_results.txt, R scripts) ---
if [ "${#ngsf_jobid[@]}" -gt 0 ]; then
  dep_list=$(IFS=:; echo "${ngsf_jobid[*]}")
  sbatch --parsable --dependency=afterany:"$dep_list" \
    --job-name=ngsF_finalize --partition=long --time=04:00:00 --mem=8gb \
    --wrap="bash '${SCRIPT_DIR}/launch_ngsF.sh' --finalize-only"
  echo "Finalize job queued (afterany:${dep_list}); will run once all ngsF jobs above finish."
else
  echo "No new ngsF jobs submitted; running finalize now."
  finalize
fi
