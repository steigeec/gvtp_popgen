#!/bin/bash
# ind_het.sh
# Per-individual heterozygosity: angsd -dosaf, then folded SFS with realSFS.
# Usage: ind_het.sh <sample.bam>   (called by launch_He.sh)
BASE_DIR=/private/groups/meyerlab/esteiger

cd "${BASE_DIR}/mapping/merged/"

file1=$1
suffix=".bam"
string=${file1%"$suffix"}
echo "processing $string"

OUT="${BASE_DIR}/diversity/He/individual"
ref="${BASE_DIR}/genome/DEINVI_genome.fasta.masked"
# all autosomal sites (variable + invariable) post-paralog-filter; no MAF filter
sites="${BASE_DIR}/selection/TajD/all_sites"
# Restricted to the 23 chromosome-scale scaffolds (the natural size break in the
# assembly, consistent with 1N=12; 82% of the post-paralog-filter genome).
# angsd's -rf/-sites intersection has a large per-scaffold setup cost, so the
# full 21,073-scaffold list took ~5-6 h/sample versus ~2 min for these 23.
# Excluded scaffolds have no sites in -sites, so output is unchanged.
chrs="${BASE_DIR}/selection/TajD/top23_chrs.txt"
NTHREADS=${NTHREADS:-3}

if [ ! -s "${OUT}/${string}.saf.gz" ]; then
    angsd -i ${file1} -nThreads ${NTHREADS} -ref $ref -anc ${ref} -dosaf 1 -GL 1 \
        -out ${OUT}/${string} -sites ${sites} -rf ${chrs} \
        -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0
    echo "done with dosaf"
fi

RS="${BASE_DIR}/software/angsd/misc/realSFS"
${RS} ${OUT}/${string}.saf.idx -fold 1 -P ${NTHREADS} > ${OUT}/${string}.ml
