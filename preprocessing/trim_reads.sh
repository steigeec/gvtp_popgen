#!/bin/bash
# Adapter- and quality-trim one sample's paired reads with Trimmomatic (PE mode).
# Writes <sample>_R1_clean.fastq.gz and <sample>_R2_clean.fastq.gz (pairs kept)
# and <sample>_S_clean.fastq.gz (reads whose mate was dropped, R1 and R2
# unpaired outputs concatenated); these are the inputs to merge_fastq.sh,
# trim_summary.sh and GVTP_mapping.sh.
#
# Conda env: trimmomatic
#   conda create -n trimmomatic -c conda-forge -c bioconda trimmomatic=0.39
# Usage: trim_reads.sh <sample>_R1_001.fastq.gz <out_dir> [threads]
set -euo pipefail

R1=$1
OUT_DIR=$2
THREADS=${3:-4}
R2=${R1/_R1_001/_R2_001}
SAMPLE=$(basename "${R1%_R1_001.fastq.gz}")

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate trimmomatic
ADAPTERS="${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE.fa"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

trimmomatic PE -threads "${THREADS}" -phred33 \
    "${R1}" "${R2}" \
    "${SAMPLE}_R1_clean.fastq.gz" "${SAMPLE}_1U.fastq.gz" \
    "${SAMPLE}_R2_clean.fastq.gz" "${SAMPLE}_2U.fastq.gz" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 LEADING:3 TRAILING:21 MINLEN:35 \
    2> "${SAMPLE}_trimmomatic.log"

cat "${SAMPLE}_1U.fastq.gz" "${SAMPLE}_2U.fastq.gz" > "${SAMPLE}_S_clean.fastq.gz"
rm "${SAMPLE}_1U.fastq.gz" "${SAMPLE}_2U.fastq.gz"
