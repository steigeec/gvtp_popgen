#!/bin/bash
# Run angsd_w_invar.sh on every genome chunk (variant + invariant sites, for Tajima's D).

#SBATCH --job-name=angsd_invar
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --mem=1000gb
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=angsd_invar_%j.log

export BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/variant_calling/angsd/chromo_30_chunks

find . -name "input_split_*" | parallel -j 1 '
  outfile=${BASE_DIR}/selection/TajD/{/.}.log
  MAF=${BASE_DIR}/selection/TajD/{/.}.mafs.gz
  if [ ! -f "$MAF" ]; then
    BAMS="${BASE_DIR}/mapping/merged/bamlist.txt"
    genome="${BASE_DIR}/genome/DEINVI_genome.fasta.masked"
    ${BASE_DIR}/gvtp_popgen/selection/TajD/angsd_w_invar.sh {} $BAMS $genome &> "$outfile"
  fi
'
