#!/bin/bash

ID=$1

cd /public/groups/meyerlab/esteiger/PreProcessing/trimmed_reads || exit 1

# Find files matching the patterns
filesR1=$(find . -type f -name "${ID}*_R1_clean.fastq.gz")
filesR2=$(find . -type f -name "${ID}*_R2_clean.fastq.gz")
filesS=$(find . -type f -name "${ID}*_S_clean.fastq.gz")

echo -e "Found R1 files: ${filesR1}"
echo -e "Found R2 files: ${filesR2}"
echo -e "Found S files: ${filesS}"

# Concatenate files
cat $filesR1 >> "${ID}_R1_clean.fastq.gz"
cat $filesR2 >> "${ID}_R2_clean.fastq.gz"
cat $filesS >> "${ID}_S_clean.fastq.gz"

# Verify the concatenated files and remove originals if valid
if [[ -f "${ID}_R1_clean.fastq.gz" && -f "${ID}_R2_clean.fastq.gz" && -f "${ID}_S_clean.fastq.gz" ]] && \
   gzip -t "${ID}_R1_clean.fastq.gz" 2>/dev/null && \
   gzip -t "${ID}_R2_clean.fastq.gz" 2>/dev/null && \
   gzip -t "${ID}_S_clean.fastq.gz" 2>/dev/null; then
    echo "All clean files for $ID are valid. Removing original files..."

    # Iterate over each list of files and remove them
    for file in $filesR1; do
        rm -f "$file" && echo "Removed $file"
    done
    for file in $filesR2; do
        rm -f "$file" && echo "Removed $file"
    done
    for file in $filesS; do
        rm -f "$file" && echo "Removed $file"
    done
else
    echo "Clean files for $ID are missing or invalid. Skipping removal."
fi
