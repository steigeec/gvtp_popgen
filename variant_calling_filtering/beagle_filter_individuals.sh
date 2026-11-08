#!/bin/bash
# Subset a beagle GL file (optionally gzipped) to the individuals listed in a
# file; the first 3 columns are always kept.
# Usage: beagle_filter_individuals.sh <inds.txt> <in.beagle[.gz]> <out.beagle>

INDS=$1
BEAGLE1=$2
BEAGLE2=$3

echo $INDS

if [[ "$BEAGLE1" == *.gz ]]; then
    read_beagle() { zcat "$BEAGLE1"; }
else
    read_beagle() { cat "$BEAGLE1"; }
fi

# Extract the header and identify matching columns
header=$(read_beagle | head -n 1)
cols=$(echo "$header" | tr '\t' '\n' | awk -v inds="$INDS" '
    BEGIN {
        while ((getline line < inds) > 0) {
            keep[line] = 1
        }
    }
    NR <= 3 {print NR}  # marker, allele1, allele2
    NR > 3 && $0 in keep {print NR}
')

read_beagle | awk -v cols="$cols" '
    BEGIN {
        split(cols, c, "\n")
        for (i in c) keep[c[i]] = 1
    }
    {
        first = 1
        for (i = 1; i <= NF; i++) {
            if (i in keep) {
                if (first) {
                    printf "%s", $i
                    first = 0
                } else {
                    printf "\t%s", $i
                }
            }
        }
        printf "\n"
    }
' > "$BEAGLE2"
