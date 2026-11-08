#!/bin/bash

#SBATCH --job-name=angsd_TajD_snps
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --mem=200gb
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=angsd_TajD_snps_%j.log

# Run variant_call_noMAF.sh (ANGSD) on each genome chunk without existing output.

export BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/variant_calling/angsd/chromo_30_chunks

find . -name "input_split_*" | parallel -j 2 '
  outfile=${BASE_DIR}/variant_calling/angsd/{/.}.log
  MAF=${BASE_DIR}/variant_calling/angsd/{/.}.mafs.gz
  if [ ! -f "$MAF" ]; then
    BAMS="${BASE_DIR}/mapping/merged/bamlist.txt"
    genome="${BASE_DIR}/genome/DEINVI_genome.fasta.masked"
    ${BASE_DIR}/gvtp_popgen/variant_calling_filtering/variant_call_noMAF.sh {} $BAMS $genome &> "$outfile"
  fi
'
