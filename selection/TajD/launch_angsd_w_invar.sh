#!/bin/bash

#SBATCH --job-name=angsd_invar
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --mem=1000gb
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=angsd_invar_%j.log

cd /private/groups/meyerlab/esteiger/variant_calling/angsd/chromo_30_chunks

find . -name "input_split_*" | parallel -j 1 '
  outfile=/private/groups/meyerlab/esteiger/selection/TajD/{/.}.log
  MAF=/private/groups/meyerlab/esteiger/selection/TajD/{/.}.mafs.gz
  if [ ! -f "$MAF" ]; then
    BAMS="/private/groups/meyerlab/esteiger/mapping/merged/bamlist.txt"
    genome="/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked"
    /private/groups/meyerlab/esteiger/selection/TajD/angsd_w_invar.sh {} $BAMS $genome &> "$outfile"
  fi
'
