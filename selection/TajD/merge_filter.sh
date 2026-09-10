#!/bin/bash
#SBATCH --job-name=TajD_merge
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=esteiger@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --output=TajD_merge_%j.log

TAJD_DIR="/private/groups/meyerlab/esteiger/selection/TajD"
CHUNKS_DIR="/private/groups/meyerlab/esteiger/variant_calling/angsd/chromo_30_chunks"

cd "$TAJD_DIR"

# ---------------------------------------------------------------------------
# COMPLETED: Pre-check and merge ran successfully on 2026-03-28.
# all_sites.mafs.gz and all_sites_filter.mafs (unfiltered copy) already exist.
#
# echo "=== Pre-check: verifying all chunks completed successfully ==="
# missing=0
# for chunk in "$CHUNKS_DIR"/input_split_*; do
#   name=$(basename "$chunk")
#   arg="${TAJD_DIR}/${name}.arg"
#   if [ ! -f "$arg" ]; then
#     echo "  MISSING .arg: ${name}"
#     missing=$((missing + 1))
#   elif ! tail -1 "$arg" | grep -q "\[ALL done\] walltime used ="; then
#     echo "  INCOMPLETE: ${name}"
#     missing=$((missing + 1))
#   fi
# done
# if [ "$missing" -gt 0 ]; then
#   echo "ERROR: $missing chunk(s) did not complete successfully. Aborting."
#   exit 1
# fi
# echo "  All chunks completed successfully. Proceeding."
#
# echo "=== Merging per-chunk mafs.gz into all_sites.mafs ==="
# header_file=$(ls input_split_*.mafs.gz | sort | head -1)
# zcat "$header_file" | head -1 > tmp_header.txt
# find . -name "input_split_*.mafs.gz" | sort | xargs zcat \
#   | tail -n +2 \
#   | sort -k1,1 -k2,2n | uniq > tmp_data.txt
# cat tmp_header.txt tmp_data.txt > all_sites.mafs
# rm tmp_header.txt tmp_data.txt
# echo "  -> $(wc -l < all_sites.mafs) lines (including header) in all_sites.mafs"
#
# echo "=== Running VerifyPolyV11 ==="
# cp all_sites.mafs all_sites_filter.mafs  # all_sites_filter.mafs already exists
# ---------------------------------------------------------------------------

#echo "=== Running VerifyPolyV11 ==="
#./VerifyPolyV11.py .999 all_sites_filter

# ---------------------------------------------------------------------------
echo ""
echo "=== Summary ==="
before=$(( $(wc -l < all_sites.mafs) - 1 ))
after=$(( $(wc -l < all_sites_filter.mafs) - 1 ))
paralog=$(wc -l < /private/groups/meyerlab/esteiger/variant_calling/angsd/postParalog.sites)
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
/private/groups/meyerlab/esteiger/software/angsd_dev/angsd sites index all_sites_filter
echo "  Created and indexed all_sites_filter"
echo "  Chromosomes in all_chrs_filter.txt: $(wc -l < all_chrs_filter.txt)"

# ---------------------------------------------------------------------------
# Now, create SAF file

cd /private/groups/meyerlab/esteiger/mapping/merged
outdir=/private/groups/meyerlab/esteiger/selection/TajD
sites_file=/private/groups/meyerlab/esteiger/selection/TajD/all_sites_filter
regions=/private/groups/meyerlab/esteiger/selection/TajD/all_chrs_filter.txt
REF=/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked
# I only want bams from a SINGLE POP 
BAMLIST="/private/groups/meyerlab/esteiger/mapping/merged/bamlist_4x_strauss.txt"
ANGSD="/private/groups/meyerlab/esteiger/software/angsd_dev/angsd"
 
$ANGSD -nThreads 15 -bam ${BAMLIST} -ref $REF -anc $REF -out $outdir/4x_strauss -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 20 -minQ 20 -minInd 50 -setMinDepthInd 1 -setMaxDepthInd 25 -doSaf 1 -rf $regions -sites $sites_file -GL 2
 
cd $outdir
RS=/private/groups/meyerlab/esteiger/software/angsd_dev/misc/realSFS
$RS 4x_strauss.saf.idx -fold 1 -P 15 > 4x_strauss.ml
 
# Calculate the thetas for each site
$RS saf2theta 4x_strauss.saf.idx -sfs 4x_strauss.ml -outname TajD
# The output from the above command are two files out.thetas.gz and out.thetas.idx.
# formal description can be found in the doc/formats.pdf in the angsd package.
 
# Estimate Tajimas D and other statistics
# win/step chosen with Rasmus
#/private/groups/meyerlab/esteiger/software/angsd_dev/misc/thetaStat do_stat TajD.thetas.idx -win 400 -step 100 -outnames TajD.thetasWindow.400
 
# Instead of using ANGSD sliding windows, where window size is by count of all bps
# try using Diana's script, which only counts SNPs
/private/groups/meyerlab/esteiger/software/angsd_dev/misc/thetaStat print TajD.thetas.idx > TajD.thetas.persite.txt
 
# Usage  	Tajimas_SNPmidpoint.py inFile outFile nSites nIndiv
# inFile 	name of infile, thetaStat print output( $pop.thetas.persite.txt )
# outFile	name of output file
# nSites 	number of sites per window
# nIndiv 	number of individuals
./Tajimas_SNPmidpoint.py TajD.thetas.persite.txt TajD.Diana.out 100 70



