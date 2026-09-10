#!/bin/bash

BAM=$1
OUTPUT_FILE=$2

CP_LIST="/private/groups/meyerlab/esteiger/genome/chloroplast.txt"
MITO_LIST="/private/groups/meyerlab/esteiger/genome/mito.txt"

prefix=$(basename "$BAM" | cut -d'_' -f1)

# Get per-scaffold sizes from BAM header, compute category sizes
samtools view -H "$BAM" | grep '^@SQ' | \
    awk -F'\t' '{
        split($2, a, ":"); split($3, b, ":");
        print a[2], b[2]
    }' > /tmp/${prefix}_scaffolds.txt

CP_SIZE=$(grep -Ff "$CP_LIST" /tmp/${prefix}_scaffolds.txt | awk '{sum+=$2} END {print sum+0}')
MITO_SIZE=$(grep -Ff "$MITO_LIST" /tmp/${prefix}_scaffolds.txt | awk '{sum+=$2} END {print sum+0}')
TOTAL_SIZE=$(awk '{sum+=$2} END {print sum+0}' /tmp/${prefix}_scaffolds.txt)
AUTO_SIZE=$((TOTAL_SIZE - CP_SIZE - MITO_SIZE))

rm /tmp/${prefix}_scaffolds.txt

# Single-pass depth: classify each position as cp, mito, or auto
read CP_SUM MITO_SUM AUTO_SUM < <(
    samtools depth "$BAM" | awk \
        -v cp_list="$CP_LIST" \
        -v mito_list="$MITO_LIST" \
    'BEGIN {
        while ((getline line < cp_list) > 0)   cp[line]   = 1
        while ((getline line < mito_list) > 0) mito[line] = 1
    }
    {
        if      ($1 in cp)   cp_sum   += $3
        else if ($1 in mito) mito_sum += $3
        else                 auto_sum += $3
    }
    END { print cp_sum+0, mito_sum+0, auto_sum+0 }'
)

TOTAL_SUM=$((AUTO_SUM + CP_SUM + MITO_SUM))

TOTAL_CVG=$(bc <<< "scale=4; $TOTAL_SUM / $TOTAL_SIZE")
AUTO_CVG=$(bc  <<< "scale=4; $AUTO_SUM  / $AUTO_SIZE")
CP_CVG=$(bc    <<< "scale=4; $CP_SUM    / $CP_SIZE")
MITO_CVG=$(bc  <<< "scale=4; $MITO_SUM  / $MITO_SIZE")

echo -e "$prefix\t$TOTAL_CVG\t$AUTO_CVG\t$CP_CVG\t$MITO_CVG" >> "$OUTPUT_FILE"
