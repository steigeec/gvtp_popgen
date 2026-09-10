#!/bin/bash

INDS=$1
BEAGLE1=$2
BEAGLE2=$3

echo $INDS

# Support gzip-compressed input
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
            keep[line] = 1  # Store individuals from $INDS
        }
    }
    NR <= 3 {print NR}  # Always keep the first 3 columns
    NR > 3 && $0 in keep {print NR}  # Add matching columns
')

# Filter and create $BEAGLE2
read_beagle | awk -v cols="$cols" '
    BEGIN {
        split(cols, c, "\n")  # Parse column indices into an array
        for (i in c) keep[c[i]] = 1
    }
    {
        first = 1
        for (i = 1; i <= NF; i++) {
            if (i in keep) {
                if (first) {
                    printf "%s", $i  # Print first column without leading tab
                    first = 0
                } else {
                    printf "\t%s", $i  # Add a tab before subsequent columns
                }
            }
        }
        printf "\n"
    }
' > "$BEAGLE2"
