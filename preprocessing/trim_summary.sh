#!/bin/bash
# Read counts before and after trimming, fraction retained, and mean trimmed
# read length for one sample; prints one tab-separated row.
# Usage: trim_summary.sh <R1_clean.fastq.gz> <dir_name>

file1=$1
dir_name=$2

# establish paths for trimmed files
file2=${file1//R1/R2}
file3=${file1//R1/S}
suffix="_R1_clean.fastq.gz"
name=${file1%"$suffix"}

# matching raw fastq files
R1raw=${file1//_R1_clean/_R1_001}
R1raw="/redser5/raw/Meyer_GVTP_reseq/${R1raw}"
R2raw=${R1raw//R1/R2}

R1raw=$(awk '{s++} END {print s/4}' "$R1raw")
R1raw=$(echo "$R1raw" | awk '{printf("%d\n", $0)}')

R2raw=$(awk '{s++} END {print s/4}' "$R2raw")
R2raw=$(echo "$R2raw" | awk '{printf("%d\n", $0)}')

R1trim=$(awk '{s++} END {print s/4}' "$file1")
R1trim=$(echo "$R1trim" | awk '{printf("%d\n", $0)}')

R2trim=$(awk '{s++} END {print s/4}' "$file2")
R2trim=$(echo "$R2trim" | awk '{printf("%d\n", $0)}')

Strim=$(awk '{s++} END {print s/4}' "$file3")
Strim=$(echo "$Strim" | awk '{printf("%d\n", $0)}')

R1perc=$(echo "scale=2; $R1trim / $R1raw" | bc)
R2perc=$(echo "scale=2; $R2trim / $R2raw"| bc)

if (( R1raw + R2raw != 0 )); then
total=$(echo "scale=2; ($R1trim + $R2trim + ($Strim*2))/($R1raw+$R2raw)" | bc)
fi

# calculate average length of reads for a sample after trimming
length=$(awk 'BEGIN {total=0; count=0}
     {if(NR%4==2) {total+=length($1); count++}}
     END {if(count>0) print total/count}' "$file1" "$file2" "$file3")


echo -e "${name}\t${R1raw}\t${R2raw}\t${R1trim}\t${R2trim}\t${Strim}\t${R1perc}\t${R2perc}\t${total}\t${length}"

