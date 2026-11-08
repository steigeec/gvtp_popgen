#!/bin/bash
# launch_He.sh
# Run ind_het.sh (per-individual SAF and folded SFS) on every BAM in bamlist.txt
# with GNU parallel. Samples with a non-empty .ml are skipped, and lock
# directories make concurrent submissions safe.
# Usage: sbatch launch_He.sh
#SBATCH --job-name=He_ind
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=90    # must equal PARALLEL_JOBS * NTHREADS
# realSFS memory scales with retained sites (up to ~330M sites, ~10 GB for one
# worker); 40gb for 30 workers caused cgroup OOM kills.
#SBATCH --mem=200gb
#SBATCH --time=168:00:00
#SBATCH --output=He_ind_%j.log

# angsd's -sites/-rf intersection sustains only ~1.5-2 cores per worker
# regardless of -nThreads, so run more workers at fewer threads each.
# Keep --cpus-per-task = PARALLEL_JOBS * NTHREADS.
BASE_DIR=/private/groups/meyerlab/esteiger

PARALLEL_JOBS=30
export NTHREADS=3

cd "${BASE_DIR}/mapping/merged/"
BAM="bamlist.txt"
OUT_DIR="${BASE_DIR}/diversity/He/individual"
SCRIPT_PATH="${BASE_DIR}/gvtp_popgen/diversity/heterozygosity/individual/ind_het.sh"

shuf ${BAM} | parallel -j ${PARALLEL_JOBS} '
  full_path={}
  filename=$(basename "$full_path")
  string=${filename%.bam}
  mlfile="'"$OUT_DIR"'"/${string}.ml
  lockdir="'"$OUT_DIR"'"/${string}.lock

  # -s so a 0-byte .ml from a failed run is not treated as complete
  if [ -s "$mlfile" ]; then
    exit 0
  fi

  # Atomic lock, safe across concurrent submissions
  if ! mkdir "$lockdir" 2>/dev/null; then
    exit 0
  fi

  # Re-check after acquiring lock
  if [ -s "$mlfile" ]; then
    rmdir "$lockdir" 2>/dev/null
    exit 0
  fi

  logfile="'"$OUT_DIR"'"/${string}_run.log
  echo "Processing: $string"
  "'"$SCRIPT_PATH"'" "$full_path" &> "$logfile"
  rmdir "$lockdir" 2>/dev/null
'
