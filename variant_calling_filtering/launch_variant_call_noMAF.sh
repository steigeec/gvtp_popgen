#!/bin/bash

#SBATCH --job-name=angsd_TajD_snps
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=esteiger@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --mem=200gb
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=angsd_TajD_snps_%j.log

cd /private/groups/meyerlab/esteiger/variant_calling/angsd/chromo_30_chunks

find . -name "input_split_*" | parallel -j 2 '
  outfile=/private/groups/meyerlab/esteiger/variant_calling/angsd/{/.}.log
  MAF=/private/groups/meyerlab/esteiger/variant_calling/angsd/{/.}.mafs.gz
  if [ ! -f "$MAF" ]; then
    BAMS="/private/groups/meyerlab/esteiger/mapping/merged/bamlist.txt"
    genome="/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked"
    /private/groups/meyerlab/esteiger/variant_calling/angsd/variant_call_noMAF.sh {} $BAMS $genome &> "$outfile"
  fi
'
