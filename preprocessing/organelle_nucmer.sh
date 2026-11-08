#!/bin/bash
# Align every genome scaffold to published Asteraceae chloroplast and
# mitochondrial genomes with NUCmer, to identify organellar scaffolds
# (classified by organelle_scaffolds.py).
# References: Helianthus annuus mitogenome KF815390.1, plastome NC_007977.1;
#   Arctium tomentosum mitogenome NC_058643.1, plastome MZ416781.1
#   (downloaded from NCBI nuccore into BASE_DIR/genome/organelle_refs/).
# Output: BASE_DIR/genome/organelle_refs/organelle_vs_genome.coords
# Conda env: mummer
#   conda create -n mummer -c conda-forge -c bioconda mummer4=4.0.0 samtools
# Usage: sbatch organelle_nucmer.sh
#SBATCH --job-name=organelle_nucmer
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/genome/organelle_refs/organelle_nucmer_%j.log
set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

REF_DIR=${BASE_DIR}/genome/organelle_refs
GENOME=${BASE_DIR}/genome/DEINVI_genome.fasta.masked
MUMMER=/private/home/esteiger/miniconda3/envs/mummer/bin

cd "${REF_DIR}"

# One reference file; prefix each record with its organelle for the classifier
for acc in NC_007977.1 MZ416781.1; do
    [ -s "${acc}.fa" ] || curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" -o "${acc}.fa"
    sed "s/^>/>plastid|/" "${acc}.fa"
done > organelle_refs.fa
for acc in KF815390.1 NC_058643.1; do
    [ -s "${acc}.fa" ] || curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" -o "${acc}.fa"
    sed "s/^>/>mito|/" "${acc}.fa"
done >> organelle_refs.fa

# Organelle genomes as reference, scaffolds as query
"${MUMMER}/nucmer" -t "${SLURM_CPUS_PER_TASK:-16}" -p organelle_vs_genome organelle_refs.fa "${GENOME}"
"${MUMMER}/show-coords" -rclTH organelle_vs_genome.delta > organelle_vs_genome.coords
echo "Alignments: $(wc -l < organelle_vs_genome.coords)"
