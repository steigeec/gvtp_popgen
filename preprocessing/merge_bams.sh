#!/bin/bash
# Merge all mapped BAMs for one sample, set read groups, remove duplicates
# (Picard MarkDuplicates) and index; writes ../merged/<ID>_rmdup.bam.
# Usage: merge_bams.sh <sample_ID>  (run by parallelized_merging.sh from mapping/)

BASE_DIR=/private/groups/meyerlab/esteiger

I="$1"
picard="${BASE_DIR}/software/picard.jar"
export TMPDIR="${BASE_DIR}/tmp/"

# strip surrounding whitespace
I=$(echo "$I" | xargs)

check_bam_validity() {
    local bamfile="$1"
    if [[ -f "$bamfile" && -s "$bamfile" ]]; then
        # samtools quickcheck if available, otherwise accept an indexed BAM
        if command -v samtools &> /dev/null && samtools quickcheck "$bamfile"; then
            return 0  # BAM is valid
        elif [[ -f "${bamfile}.bai" ]]; then
            return 0  # BAM is assumed complete if indexed
        fi
    fi
    return 1  # BAM is invalid or incomplete
}

# skip samples that are already done
if check_bam_validity "${BASE_DIR}/merged/${I}_rmdup.bam"; then
    exit 0
else
    echo "No valid BAM file found. Proceeding with the script."
fi

echo "Processing sample ID: '$I'"
echo "Searching for BAM files pertaining to sample $I"

declare -a bam_files=()
declare -a bai_files=()

while IFS= read -r file; do
  echo "Found BAM file: $file"
  bam_files+=("$file")
done < <(find . -type f -name "*${I}*dup.bam")

while IFS= read -r file; do
  echo "Found BAI file: $file"
  bai_files+=("$file")
done < <(find . -type f -name "*${I}*dup.bai")

bam_count=${#bam_files[@]}

if [ "$bam_count" -eq 0 ]; then
  echo "No BAM files found for sample $I"
elif [ "$bam_count" -eq 1 ]; then
  echo "Only one BAM file found for sample $I. Copying the file."
  echo "Fixing read groups"
  java -Djava.io.tmpdir="${BASE_DIR}/tmp/" -Xmx6g \
    -jar $picard AddOrReplaceReadGroups \
    I="${bam_files[0]}" \
    O="../merged/${I}.bam" \
    RGID=${I} \
    RGLB=${I} \
    RGPL=ILLUMINA \
    RGPU=${I} \
    RGSM=${I}
  if [ $? -eq 0 ]; then
    echo "Successfully copied ${bam_files[0]} to ../merged/${I}.bam"
    java -Djava.io.tmpdir="${BASE_DIR}/tmp/" -Xmx6g -jar $picard MarkDuplicates \
	REMOVE_DUPLICATES=true M=../merged/${I}.metrics \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1300 I=../merged/${I}.bam \
	O=../merged/${I}_rmdup.bam
    samtools index ../merged/${I}_rmdup.bam && \
           rm "../merged/${I}.bam" "${bam_files[@]}" && \
           find . -name "*${I}*_markdup.bai" -exec rm {} +
  else
    echo "Failed to copy BAM file for sample $I"
  fi
else
  echo "Merging BAM files for sample $I"
  samtools merge -f -@2 "../merged/${I}.bam" "${bam_files[@]}"
  echo "Fixing read groups"
java -Djava.io.tmpdir="${BASE_DIR}/tmp/" -Xmx7g \
    -jar $picard AddOrReplaceReadGroups \
    I="../merged/${I}.bam" \
    O="../merged/${I}_RGcorrected.bam" \
    RGID=${I} \
    RGLB=${I} \
    RGPL=ILLUMINA \
    RGPU=${I} \
    RGSM=${I}  
  if [ $? -eq 0 ]; then
    echo "Successfully merged BAM files into ../merged/${I}.bam"
    java -Djava.io.tmpdir="${BASE_DIR}/tmp/" -Xmx7g \
	-jar $picard MarkDuplicates REMOVE_DUPLICATES=true \
	M=../merged/${I}.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1300 \
	I="../merged/${I}_RGcorrected.bam" O="../merged/${I}_rmdup.bam"
    samtools index ../merged/${I}_rmdup.bam && \
           rm "../merged/${I}.bam" "${bam_files[@]}" "../merged/${I}_RGcorrected.bam" && \
           find . -name "*${I}*_markdup.bai" -exec rm {} +
  else
    echo "Failed to merge BAM files for sample $I"
  fi
fi