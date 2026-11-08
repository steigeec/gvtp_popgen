#!/bin/bash
# Run small_chromplot.py with no interval first, then retry with increasing -i
# values until one attempt finishes within PER_ATTEMPT_TIMEOUT.
# Usage: plot_with_retry.sh <prefix> <chrom:start-end> <out.png> <legend.png> [small_chromplot.py args]

set -u

BASE_DIR=/private/groups/meyerlab/esteiger

PREFIX=$1        # e.g. DEINVI_00006_GL_maf5
REGION=$2        # e.g. DEINVI_00006:500000-70016456
OUT=$3            # output png path
LEGEND_OUT=$4     # legend png path
shift 4
EXTRA_ARGS=("$@") # e.g. -m metadata.tsv -g fog --numeric

SCRIPT=${BASE_DIR}/gvtp_popgen/inversion/small_chromplot.py
PER_ATTEMPT_TIMEOUT=600   # 10 min per attempt

for interval in "" 5 10 15 20 25 30; do
    if [ -f "$OUT" ]; then
        echo "Output already exists, skipping remaining attempts: $OUT"
        break
    fi

    if [ -z "$interval" ]; then
        IARG=()
        TAG="no-interval"
    else
        IARG=(-i "$interval")
        TAG="interval-$interval"
    fi

    echo "Attempt: $TAG (timeout ${PER_ATTEMPT_TIMEOUT}s)"
    timeout ${PER_ATTEMPT_TIMEOUT} python3 "$SCRIPT" "$PREFIX" "$REGION" \
        "${EXTRA_ARGS[@]}" -o "$OUT" --legend_out "$LEGEND_OUT"
    STATUS=$?

    if [ -f "$OUT" ]; then
        echo "Succeeded at $TAG"
        break
    else
        echo "Failed/timed out at $TAG (exit $STATUS)"
    fi
done

if [ ! -f "$OUT" ]; then
    echo "All attempts exhausted without producing $OUT"
    exit 1
fi
