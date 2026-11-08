#!/bin/bash
# Concatenate a sample's trimmed R1/R2/S fastq files across sequencing runs and
# remove the per-run files once the merged files pass gzip -t.
# Usage: merge_fastq.sh <sample_ID>

ID=$1

cd /public/groups/meyerlab/esteiger/PreProcessing/trimmed_reads || exit 1

filesR1=$(find . -type f -name "${ID}*_R1_clean.fastq.gz")
filesR2=$(find . -type f -name "${ID}*_R2_clean.fastq.gz")
filesS=$(find . -type f -name "${ID}*_S_clean.fastq.gz")

echo -e "Found R1 files: ${filesR1}"
echo -e "Found R2 files: ${filesR2}"
echo -e "Found S files: ${filesS}"

cat $filesR1 >> "${ID}_R1_clean.fastq.gz"
cat $filesR2 >> "${ID}_R2_clean.fastq.gz"
cat $filesS >> "${ID}_S_clean.fastq.gz"

# Verify the concatenated files and remove originals if valid
if [[ -f "${ID}_R1_clean.fastq.gz" && -f "${ID}_R2_clean.fastq.gz" && -f "${ID}_S_clean.fastq.gz" ]] && \
   gzip -t "${ID}_R1_clean.fastq.gz" 2>/dev/null && \
   gzip -t "${ID}_R2_clean.fastq.gz" 2>/dev/null && \
   gzip -t "${ID}_S_clean.fastq.gz" 2>/dev/null; then
    echo "All clean files for $ID are valid. Removing original files..."

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
