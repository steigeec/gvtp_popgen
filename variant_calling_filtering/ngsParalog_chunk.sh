#!/bin/bash

#SBATCH --job-name=ngsParalog_chunk
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --time=168:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/variant_calling/ngsParalog/logs/chunk_%A_%a.log

# ngsParalog calcLR on one chunk of sites (array task), via samtools mpileup.
# Usage: submitted by launch_ngsParalog.sh

BASE_DIR=/private/groups/meyerlab/esteiger

CHUNKS_DIR="$1"
BAMLIST="${BASE_DIR}/mapping/merged/bamlist.txt"
NGSPARALOG="${BASE_DIR}/software/ngsParalog/ngsParalog"
SAMTOOLS="${BASE_DIR}/software/samtools/bin/samtools"

mapfile -t chunks < <(find "${CHUNKS_DIR}" -name "sites_chunk_*" -not -name "*.*" | sort)
sites="${chunks[$SLURM_ARRAY_TASK_ID]}"
outfile="${sites}.lr"
lockdir="${sites}.lock"

if [ -f "$outfile" ]; then
  echo "Output already exists for ${sites}, skipping."
  exit 0
fi

# mkdir is atomic, so only one task can take the lock
if ! mkdir "$lockdir" 2>/dev/null; then
  echo "Another task is processing ${sites}, skipping."
  exit 0
fi
trap "rmdir '$lockdir'" EXIT

echo "Processing ${sites}..."
bed="${sites}.bed"
pileup="${sites}.pileup"
trap "rmdir '$lockdir' 2>/dev/null; rm -f '$bed'" EXIT

if [ -f "$pileup" ]; then
  echo "Pileup already exists, skipping samtools step."
else
  awk '{print $1 "\t" ($2-1) "\t" $2}' "$sites" > "$bed"
  EXPECTED_LINES=$(wc -l < "$sites")
  cd ${BASE_DIR}/mapping/merged
  # -minQ 20 in ngsParalog, as recommended on the ngsParalog GitHub
  # -q 0 -Q 0 for mpileup, as recommended on the ngsParalog GitHub
  # pileup goes to a file, not a pipe: with 842 samples one line can exceed the
  #   64 KB pipe buffer and samtools and ngsParalog deadlock
  # -d 110000 so samtools does not silently downsample high-depth sites
  "$SAMTOOLS" mpileup \
    --ff UNMAP,QCFAIL,DUP \
    -q 0 -Q 0 -d 110000 \
    -b "$BAMLIST" \
    -l "$bed" \
    > "$pileup" &
  SAMTOOLS_PID=$!
  # samtools can stall at the very end: if within 9 lines of expected and idle 5 h, treat as done
  while kill -0 "$SAMTOOLS_PID" 2>/dev/null; do
    sleep 900
    CURRENT_LINES=$(wc -l < "$pileup" 2>/dev/null || echo 0)
    DIFF=$((EXPECTED_LINES - CURRENT_LINES))
    if [ "$DIFF" -ge 0 ] && [ "$DIFF" -le 9 ]; then
      AGE=$(( $(date +%s) - $(stat -c %Y "$pileup" 2>/dev/null || echo 0) ))
      if [ "$AGE" -ge 18000 ]; then
        echo "Pileup within ${DIFF} lines of expected, no writes for ${AGE}s — assuming complete, killing samtools."
        kill "$SAMTOOLS_PID" 2>/dev/null
        break
      fi
    fi
  done
  wait "$SAMTOOLS_PID" 2>/dev/null
fi

"$NGSPARALOG" calcLR \
  -infile "$pileup" -minQ 20 \
  -outfile "$outfile"
ngs_exit=$?

if [ $ngs_exit -eq 0 ] && [ -s "$outfile" ]; then
  rm -f "$pileup"
else
  echo "WARNING: ngsParalog exited ${ngs_exit} or .lr missing/empty — keeping pileup for retry"
fi
rm -f "$bed"
