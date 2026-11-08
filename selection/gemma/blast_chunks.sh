#!/bin/bash
# BLAST one chunk of candidate sequences against remote NCBI nt (megablast).

#SBATCH --job-name=blast_chunk
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=12:00:00
#SBATCH --output=blast_chunk_%A_%a.log

# Usage: sbatch --array=0-N blast_chunks.sh <var_dir>
# e.g.: sbatch --array=0-6 blast_chunks.sh BASE_DIR/selection/gemma/fog
# Conda env: gemoma (BASE_DIR/software/miniconda3/envs/gemoma; for blastn)
#   conda create -n gemoma -c bioconda gemoma=1.9 mmseqs2
#   (BLAST+ 2.16.0 comes in as a GeMoMa dependency)

BASE_DIR=/private/groups/meyerlab/esteiger

VAR_DIR=$1
CHUNK=$(printf "chunk_%02d" ${SLURM_ARRAY_TASK_ID})

echo "BLASTing ${CHUNK} in ${VAR_DIR}"

${BASE_DIR}/software/miniconda3/envs/gemoma/bin/blastn -query "${VAR_DIR}/blast_chunks/${CHUNK}.fasta" \
       -db nt \
       -remote \
       -max_target_seqs 15 \
       -evalue 1e-10 \
       -task megablast \
       -outfmt "6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
       -out "${VAR_DIR}/blast_chunks/${CHUNK}.out"

echo "Done: ${CHUNK} — $(wc -l < "${VAR_DIR}/blast_chunks/${CHUNK}.out") hits"
