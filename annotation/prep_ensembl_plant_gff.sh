#!/bin/bash
# Download a species' whole-genome GFF3 from Ensembl Plants (release 63), for
# GENESPACE. The gff3/ directory slug can differ from the fasta/ slug (e.g.
# lactuca_sativa_gca002870075v4gb); check the FTP listing. Skips the download
# if the file already exists.
#
# Usage: bash prep_ensembl_plant_gff.sh <species_slug> <file_prefix> <display_name>
#   species_slug: Ensembl Plants FTP gff3/ directory name, e.g.
#                 lactuca_sativa_gca002870075v4gb
#   file_prefix:  the gff3 filename prefix before ".gff3.gz", e.g.
#                 "Lactuca_sativa_gca002870075v4gb.Lsat_Salinas_v11.63"
#   display_name: short name used for the output filename, e.g. "Lactuca_sativa"
#
# Example:
#   bash prep_ensembl_plant_gff.sh lactuca_sativa_gca002870075v4gb \
#       Lactuca_sativa_gca002870075v4gb.Lsat_Salinas_v11.63 Lactuca_sativa

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

if [ $# -lt 3 ]; then
    echo "Usage: bash prep_ensembl_plant_gff.sh <species_slug> <file_prefix> <display_name>" >&2
    exit 1
fi

SPECIES_SLUG="$1"
FILE_PREFIX="$2"
DISPLAY_NAME="$3"

OUTDIR=${BASE_DIR}/genome/annotation
RELEASE=63
URL="https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-${RELEASE}/gff3/${SPECIES_SLUG}/${FILE_PREFIX}.gff3.gz"

OUT_GZ="${OUTDIR}/${DISPLAY_NAME}.gff3.gz"

if [ -s "$OUT_GZ" ]; then
    echo "GFF3 for ${DISPLAY_NAME} already downloaded at ${OUT_GZ} — skipping download."
else
    echo "Downloading ${URL}..."
    curl -sS --fail -o "$OUT_GZ" "$URL"
fi

echo "Done. GFF3: ${OUT_GZ}"
