#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH --partition=long
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=20gb
#SBATCH --output=/private/groups/meyerlab/esteiger/coverage_%j.log

OUTPUT_FILE="/private/groups/meyerlab/esteiger/avg_cvg.txt"

echo -e "sample\ttotal_cvg\tauto_cvg\tcp_cvg\tmito_cvg" > "$OUTPUT_FILE"

find /private/groups/meyerlab/esteiger/mapping/merged -maxdepth 1 -name '*_rmdup.bam' | \
parallel -j10 --no-run-if-empty --bar '
    bam={}
    prefix=$(basename "$bam" | cut -d"_" -f1)
    if ! grep -q "^${prefix}[[:space:]]" "'"$OUTPUT_FILE"'"; then
        /private/groups/meyerlab/esteiger/coverage.sh "$bam" "'"$OUTPUT_FILE"'"
    fi
'
