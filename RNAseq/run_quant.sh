#!/bin/bash
# Salmon TPM for the 8 tissue RNA-seq libraries (McEvoy et al. 2024; downloaded
# from NCBI already trimmed): extract transcripts from the masked genome with the
# longest-isoform GFF, build a salmon index, quantify each library.
# Output: RNA_DIR/tissue_expression/quant/<sample>/quant.sf
# Conda env: rnaseq_quant
#   conda create -n rnaseq_quant -c conda-forge -c bioconda salmon=2.7.0 gffread=0.12.9
# Usage: bash run_quant.sh
set -euo pipefail

RNA_DIR=/scratch1/emma/RNA-seq
GENOME_DIR=/redser5/projects/GVTP/genome/gene_annotations/easel

cd "${RNA_DIR}/tissue_expression"

source /redser5/software/anaconda/etc/profile.d/conda.sh
conda activate rnaseq_quant

GENOME=${GENOME_DIR}/DEINVI_genome.masked.fasta
GFF=${GENOME_DIR}/fixes/longest_isoform.renamed.rm_internal_stops.gff
TX_FASTA=ref/transcripts.fa
INDEX=ref/salmon_index
THREADS=8

echo "=== [1/3] Extracting transcript sequences with gffread ==="
if [ ! -s "$TX_FASTA" ]; then
    gffread -w "$TX_FASTA" -g "$GENOME" "$GFF"
fi
grep -c "^>" "$TX_FASTA"

echo "=== [2/3] Building salmon index ==="
if [ ! -d "$INDEX" ]; then
    salmon index -t "$TX_FASTA" -i "$INDEX" -k 31 -p "$THREADS"
fi

echo "=== [3/3] Quantifying samples ==="
FASTQ_DIR=$(find data -type d -name trimmed | head -1)
echo "Using fastq dir: $FASTQ_DIR"

declare -a SAMPLES=(GVTP-1Root_S2 GVTP-1Shoot_S1 GVTP-2Root_S4 GVTP-2Shoot_S3 GVTP-5Root_S6 GVTP-5Shoot_S5 GVTPMLeaf_S7 GVTPMRoot_S8)

for s in "${SAMPLES[@]}"; do
    r1=$(find "$FASTQ_DIR" -name "trimmed_${s}_L002_R1_001.fastq.gz")
    r2=$(find "$FASTQ_DIR" -name "trimmed_${s}_L002_R2_001.fastq.gz")
    outdir="quant/${s}"
    if [ -s "${outdir}/quant.sf" ]; then
        echo "Skipping $s (already done)"
        continue
    fi
    echo "--- Quantifying $s ---"
    echo "R1: $r1"
    echo "R2: $r2"
    salmon quant -i "$INDEX" -l A -1 "$r1" -2 "$r2" \
        -p "$THREADS" --validateMappings -o "$outdir"
done

echo "ALL_QUANT_DONE"
