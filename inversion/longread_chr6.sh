#!/bin/bash

#SBATCH --job-name=longread_chr6
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128gb
#SBATCH --time=120:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/inversions/longreads/longread_%x_%j.log

# Align the DEINVI reference individual's long reads (SRA) to the genome and
# keep ONLY the alignments on DEINVI_00006, to check the candidate inversion
# breakpoints (e.g. ~45.3 and ~49.3 Mb for INV1).
#
# Usage: sbatch -J <name> longread_chr6.sh <SRR accession> <minimap2 preset>
#   sbatch -J hifi_chr6 longread_chr6.sh SRR27077960 map-hifi   # PacBio HiFi, SRX22767383
#   sbatch -J ont_chr6  longread_chr6.sh SRR27077961 map-ont    # ONT PromethION, SRX22767382
#
# Space: the .sra (45 / 89 GB) and FASTQ (~150-200 GB) live ONLY on the
# compute node's local disk (/data/tmp, ~14 TB) and are deleted on exit.
# Nothing is written to $HOME -- the SRA toolkit's default download location
# (~/ncbi) would overflow the home quota -- so prefetch gets an explicit
# output dir and runs from the local dir. Kept: <SRR>_DEINVI_00006.bam(.bai)
# in OUT_DIR (a few GB).

set -euo pipefail

SRR=$1
PRESET=$2
CHROM=DEINVI_00006
GENOME=/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked
OUT_DIR=/private/groups/meyerlab/esteiger/inversions/longreads
LOCAL=/data/tmp/longread_${SLURM_JOB_ID}
THREADS=${SLURM_NTASKS:-32}

mkdir -p ${LOCAL} ${OUT_DIR}
trap 'rm -rf ${LOCAL}' EXIT
cd ${LOCAL}

avail_gb=$(df -BG --output=avail ${LOCAL} | tail -1 | tr -dc '0-9')
echo "$(date): ${SRR} (${PRESET}) on $(hostname), ${avail_gb} GB free in ${LOCAL}"
if (( avail_gb < 600 )); then
    echo "ERROR: need >= 600 GB of local disk" >&2
    exit 1
fi

# 1. download (default --max-size is 20 GB, far below these runs; value in KB = 500 GB)
prefetch --max-size 524288000 -O ${LOCAL} ${SRR}
SRA_FILE=$(find ${LOCAL}/${SRR} -name "${SRR}.sra*" | head -1)
echo "$(date): downloaded $(du -h ${SRA_FILE} | cut -f1)"

# 2. SRA -> FASTQ, temp files on local disk
fasterq-dump --threads 8 --temp ${LOCAL}/tmp --outdir ${LOCAL} ${SRA_FILE}
rm -f ${SRA_FILE}
FASTQ=$(ls ${LOCAL}/${SRR}*.fastq | head -1)
echo "$(date): FASTQ $(du -h ${FASTQ} | cut -f1)"

# 3. align to the whole genome (so reads from elsewhere don't pile onto
#    DEINVI_00006), keep header + every alignment on DEINVI_00006 (primary,
#    supplementary and secondary -- split parts are what mark breakpoints).
#    -Y: soft-clip supplementary alignments so their full sequence is kept.
minimap2 -t ${THREADS} -ax ${PRESET} -Y --MD ${GENOME} ${FASTQ} \
    | awk -v c=${CHROM} '/^@/ || $3 == c' \
    | samtools sort -@ 8 -m 2G -T ${LOCAL}/sort -o ${OUT_DIR}/${SRR}_${CHROM}.bam -
samtools index ${OUT_DIR}/${SRR}_${CHROM}.bam
echo "$(date): kept $(samtools view -c ${OUT_DIR}/${SRR}_${CHROM}.bam) alignments on ${CHROM} -> ${OUT_DIR}/${SRR}_${CHROM}.bam"
