#!/bin/bash
# Summarize VerifyPolyV11.py site filtering, build and index the filtered sites file,
# then compute Tajima's D for Strauss in 100-SNP windows.
# Inputs (in $TAJD_DIR): all_sites.mafs (merged per-chunk ANGSD mafs) and
# all_sites_filter.mafs (output of ./VerifyPolyV11.py .999 all_sites_filter).
#SBATCH --job-name=TajD_merge
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --output=TajD_merge_%j.log

BASE_DIR=/private/groups/meyerlab/esteiger

TAJD_DIR="${BASE_DIR}/selection/TajD"
CHUNKS_DIR="${BASE_DIR}/variant_calling/angsd/chromo_30_chunks"

cd "$TAJD_DIR"

# ---------------------------------------------------------------------------
echo ""
echo "=== Summary ==="
before=$(( $(wc -l < all_sites.mafs) - 1 ))
after=$(( $(wc -l < all_sites_filter.mafs) - 1 ))
paralog=$(wc -l < ${BASE_DIR}/variant_calling/angsd/postParalog.sites)
echo "  Sites before VerifyPoly filter: $before"
echo "  Sites after VerifyPoly filter:  $after"
echo "  Sites in postParalog.sites:     $paralog"

# ---------------------------------------------------------------------------
echo ""
echo "=== Post-processing ==="
pigz -p 100 all_sites.mafs
echo "  Compressed all_sites.mafs"

cat all_sites_filter.mafs | cut -f1,2 | tail -n +2 > all_sites_filter
cut -f1 all_sites_filter | sort | uniq > all_chrs_filter.txt
${BASE_DIR}/software/angsd_dev/angsd sites index all_sites_filter
echo "  Created and indexed all_sites_filter"
echo "  Chromosomes in all_chrs_filter.txt: $(wc -l < all_chrs_filter.txt)"

# ---------------------------------------------------------------------------
# Site allele frequency (SAF), folded SFS and per-site thetas for Strauss

cd ${BASE_DIR}/mapping/merged
outdir=${BASE_DIR}/selection/TajD
sites_file=${BASE_DIR}/selection/TajD/all_sites_filter
regions=${BASE_DIR}/selection/TajD/all_chrs_filter.txt
REF=${BASE_DIR}/genome/DEINVI_genome.fasta.masked
# Bams from a single population only
BAMLIST="${BASE_DIR}/mapping/merged/bamlist_4x_strauss.txt"
ANGSD="${BASE_DIR}/software/angsd_dev/angsd"
 
$ANGSD -nThreads 15 -bam ${BAMLIST} -ref $REF -anc $REF -out $outdir/4x_strauss -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 50 -setMinDepthInd 1 -setMaxDepthInd 25 -doSaf 1 -rf $regions -sites $sites_file -GL 2
 
cd $outdir
RS=${BASE_DIR}/software/angsd_dev/misc/realSFS
$RS 4x_strauss.saf.idx -fold 1 -P 15 > 4x_strauss.ml
 
# Calculate the thetas for each site
$RS saf2theta 4x_strauss.saf.idx -sfs 4x_strauss.ml -outname TajD
# Writes TajD.thetas.gz and TajD.thetas.idx (format described in doc/formats.pdf in the angsd package).
 
# Optional alternative: ANGSD windows by bp (win/step chosen with Rasmus)
#BASE_DIR/software/angsd_dev/misc/thetaStat do_stat TajD.thetas.idx -win 400 -step 100 -outnames TajD.thetasWindow.400
 
# Tajima's D in windows defined by SNP count rather than by bp (Diana Aguilar's script)
${BASE_DIR}/software/angsd_dev/misc/thetaStat print TajD.thetas.idx > TajD.thetas.persite.txt
 
# Usage  	Tajimas_SNPmidpoint.py inFile outFile nSites nIndiv
# inFile 	name of infile, thetaStat print output( $pop.thetas.persite.txt )
# outFile	name of output file
# nSites 	number of sites per window
# nIndiv 	number of individuals
python3 ${BASE_DIR}/gvtp_popgen/selection/TajD/Tajimas_SNPmidpoint.py TajD.thetas.persite.txt TajD.Diana.out 100 70

