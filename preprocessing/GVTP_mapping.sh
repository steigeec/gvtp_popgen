#!/bin/bash
# Map one sample's trimmed paired and orphan reads with bwa mem, merge, sort
# and index. Skips samples that already have a valid BAM.
# Usage: GVTP_mapping.sh <R1_clean.fastq.gz> <outdir> <threads>  (run by parallelized_mapping.sh)

BASE_DIR=/private/groups/meyerlab/esteiger

file1=$1
outdir=$2
threads=$3

prefix="./"
file1=${file1#"$prefix"}
suffix="_R1_clean.fastq.gz"
root=${file1%"$suffix"}
file2=${root}_R2_clean.fastq.gz
string=$(echo "$root" | sed -e 's/SampleA_//g' -e 's/SampleB_//g')
sample="${string%%_*}"
file3=${root}_S_clean.fastq.gz

export TMPDIR="${BASE_DIR}/tmp/"

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
if check_bam_validity "${outdir}/${string}_markdup.bam"; then
    echo "Valid BAM file found. Skipping the rest of the script."
    exit 0
else
    echo "No valid BAM file found. Proceeding with the script."
fi

genome=${BASE_DIR}/genome/DEINVI_genome.fasta
RUN=$(zcat "$file1" | head -n 1 | awk -F ':' '{print $2}')
bwa=${BASE_DIR}/software/bwa/bwa

echo -e "\n\n align with bwa mem for ${string}" 
# -P: in PE mode, rescue mates by Smith-Waterman only, without forcing a proper pair
# -M: mark shorter split hits as secondary (Picard compatibility)
# one library per sample, so ID = LB = SM
# pairs
$bwa mem -t ${threads} -PM -R "@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:ILLUMINA" ${genome} ${file1} ${file2} > ${outdir}/${string}_temp_p.sam
#orphans
$bwa mem -t ${threads} -M -R "@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:ILLUMINA" ${genome} ${file3} > ${outdir}/${string}_temp_u.sam

echo -e "\n\n generate bam from sam for ${string}" 
samtools view -@ ${threads} -b -o ${outdir}/${string}_temp_p.bam ${outdir}/${string}_temp_p.sam
samtools view -@ ${threads} -b -o ${outdir}/${string}_temp_u.bam ${outdir}/${string}_temp_u.sam
echo -e "\n\n merge pairs and orphans for ${string}" 
samtools merge -c ${outdir}/${string}_temp.bam ${outdir}/${string}_temp_p.bam ${outdir}/${string}_temp_u.bam
echo -e "\n\n sort by genomic requirements for ${string}" 
samtools sort -@ ${threads} -o ${outdir}/${string}_markdup.bam -T ${TMPDIR}/${string}_ ${outdir}/${string}_temp.bam

echo -e "\n\n clean up extra files for ${string}" 
rm ${outdir}/${string}_temp_p.sam ${outdir}/${string}_temp_u.sam ${outdir}/${string}_temp_p.bam ${outdir}/${string}_temp_u.bam ${outdir}/${string}_temp.bam

echo -e "\n\nindex bam for ${string}" 
samtools index ${outdir}/${string}_markdup.bam
