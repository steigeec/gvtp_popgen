#!/bin/bash

cd "/private/groups/meyerlab/esteiger/mapping/merged/"

file1=$1
suffix=".bam"
string=${file1%"$suffix"}
echo "processing $string"

OUT="/private/groups/meyerlab/esteiger/diversity/He/individual"
ref="/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked"
# all autosomal sites (variable + invariable) post-paralog-filter; no MAF filter
sites="/private/groups/meyerlab/esteiger/selection/TajD/all_sites"
# Restricted to the 23 chromosome-scale scaffolds (the natural size break in the
# assembly, consistent with 1N=12 for this diploid species; captures 82% of the
# post-paralog-filtered genome). angsd's -rf/-sites intersection pays a large
# fixed per-scaffold setup cost, so using the full 21,073-scaffold list (mostly
# tiny unplaced fragments) made -dosaf take ~5-6 hours/sample; restricting to
# these 23 scaffolds cut that to ~2 minutes in testing, with no change to what
# -sites allows to be emitted (excluded scaffolds were never in -sites anyway).
chrs="/private/groups/meyerlab/esteiger/selection/TajD/top23_chrs.txt"
NTHREADS=${NTHREADS:-3}

if [ ! -s "${OUT}/${string}.saf.gz" ]; then
    angsd -i ${file1} -nThreads ${NTHREADS} -ref $ref -anc ${ref} -dosaf 1 -GL 1 \
        -out ${OUT}/${string} -sites ${sites} -rf ${chrs} \
        -minMapQ 20 -minQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0
    echo "done with dosaf"
fi

RS="/private/groups/meyerlab/esteiger/software/angsd/misc/realSFS"
${RS} ${OUT}/${string}.saf.idx -fold 1 -P ${NTHREADS} > ${OUT}/${string}.ml
