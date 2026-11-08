#!/bin/bash
# Build the GENESPACE working directory (genespace/peptide/, genespace/bed/)
# from the longest-isoform peptide fastas and gff annotations in
# genome/annotation/ (Deinandra: genome/).
#   wd/peptide/<genomeID>.fa   protein fasta, header == bed id
#   wd/bed/<genomeID>.bed      chr, start (0-based), end, id; tab-delimited, no header
# The id_regex per species picks the gff attribute that matches its peptide
# headers (see build_genespace_bed.py); genomeIDs match orthofinder_input/.
# Safe to re-run.
#
# Usage: bash setup_genespace_inputs.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

GENOME_DIR=${BASE_DIR}/genome
ANNOT_DIR="${GENOME_DIR}/annotation"
WD=${BASE_DIR}/inversions/synteny/genespace
BUILD_BED="${BASE_DIR}/gvtp_popgen/inversion/build_genespace_bed.py"

mkdir -p "${WD}/peptide" "${WD}/bed"

# genomeID | peptide_fasta | gff | id_regex
# (pipe-delimited because the id_regex fields contain colons)
SPECIES=(
  "Deinandra|${GENOME_DIR}/remove_duplicates_rm_internalstops_longest_isoform.pep|${GENOME_DIR}/remove_duplicates_rm_internalstops_longest_isoform.gff|mRNA_id=([^;]+)"
  "Arabidopsis_thaliana|${ANNOT_DIR}/Arabidopsis_thaliana.TAIR10.pep.longest.fa|${ANNOT_DIR}/Arabidopsis_thaliana.gff3.gz|ID=gene:([^;]+)"
  "Helianthus_annuus|${ANNOT_DIR}/Helianthus_annuus.pep.longest.fa|${ANNOT_DIR}/Helianthus_annuus.gff3.gz|ID=gene:([^;]+)"
  "Helianthus_anomalus|${ANNOT_DIR}/Helianthus_anomalus.pep.longest.fa|${ANNOT_DIR}/Helianthus_anomalus.gff.gz|locus_tag=([^;]+)"
  "Lactuca_sativa|${ANNOT_DIR}/Lactuca_sativa.pep.longest.fa|${ANNOT_DIR}/Lactuca_sativa.gff3.gz|ID=gene:([^;]+)"
  "Solanum_lycopersicum|${ANNOT_DIR}/Solanum_lycopersicum.pep.longest.fa|${ANNOT_DIR}/Solanum_lycopersicum.gff3.gz|ID=gene:([^;]+)"
)

for entry in "${SPECIES[@]}"; do
  IFS="|" read -r genome_id pep gff id_regex <<< "$entry"
  echo "=== ${genome_id} ==="

  ln -sf "${pep}" "${WD}/peptide/${genome_id}.fa"

  python3 "${BUILD_BED}" "${gff}" "${id_regex}" "${pep}" "${WD}/bed/${genome_id}.bed"
  echo
done

echo "GENESPACE input directory ready: ${WD}"
echo "  peptide/: $(ls "${WD}/peptide" | wc -l) files"
echo "  bed/:     $(ls "${WD}/bed" | wc -l) files"
