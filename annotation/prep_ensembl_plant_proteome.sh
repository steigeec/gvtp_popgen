#!/bin/bash
# Download an Ensembl Plants proteome (release 63) and keep the longest
# protein per gene (grouped by the header's gene: field), for OrthoFinder.
# Skips the download if the raw file already exists.
#
# Usage: bash prep_ensembl_plant_proteome.sh <species_slug> <file_prefix> <display_name>
#   species_slug: Ensembl Plants FTP directory name, e.g. helianthus_annuus
#   file_prefix:  the FASTA filename prefix before ".pep.all.fa.gz", e.g.
#                 "Helianthus_annuus.HanXRQr2.0-SUNRISE"
#   display_name: short name used for the output filename, e.g. "Helianthus_annuus"
#
# Example:
#   bash prep_ensembl_plant_proteome.sh helianthus_annuus \
#       Helianthus_annuus.HanXRQr2.0-SUNRISE Helianthus_annuus

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

if [ $# -lt 3 ]; then
    echo "Usage: bash prep_ensembl_plant_proteome.sh <species_slug> <file_prefix> <display_name>" >&2
    exit 1
fi

SPECIES_SLUG="$1"
FILE_PREFIX="$2"
DISPLAY_NAME="$3"

OUTDIR=${BASE_DIR}/genome/annotation
RELEASE=63
URL="https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-${RELEASE}/fasta/${SPECIES_SLUG}/pep/${FILE_PREFIX}.pep.all.fa.gz"

RAW_GZ="${OUTDIR}/${DISPLAY_NAME}.pep.all.fa.gz"
RAW="${OUTDIR}/${DISPLAY_NAME}.pep.all.fa"
LONGEST="${OUTDIR}/${DISPLAY_NAME}.pep.longest.fa"

if [ -s "$RAW_GZ" ]; then
    echo "Raw ${DISPLAY_NAME} proteome already downloaded at ${RAW_GZ} — skipping download."
else
    echo "Downloading ${URL}..."
    curl -sS --fail -o "$RAW_GZ" "$URL"
fi

gunzip -k -f "$RAW_GZ"

# Headers: ">TXID pep ... gene:GENEID ..."; keep the longest protein per gene
python3 - "$RAW" "$LONGEST" <<'PY'
import sys
import re

raw_file, out_file = sys.argv[1], sys.argv[2]

best = {}  # gene_id -> (length, tx_id, seq)
gene_id = None
tx_id = None
seq_chunks = []

def flush():
    if gene_id is None:
        return
    seq = "".join(seq_chunks)
    prev = best.get(gene_id)
    if prev is None or len(seq) > prev[0]:
        best[gene_id] = (len(seq), tx_id, seq)

with open(raw_file) as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            tx_id = line[1:].split()[0]
            m = re.search(r"gene:(\S+)", line)
            gene_id = m.group(1) if m else tx_id.split(".")[0]
            seq_chunks = []
        else:
            seq_chunks.append(line)
    flush()

with open(out_file, "w") as out:
    for gene_id in sorted(best):
        _, tx_id, seq = best[gene_id]
        out.write(f">{gene_id}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i + 60] + "\n")

print(f"Wrote {len(best)} longest-isoform gene proteins to {out_file}")
PY

echo "Done. Longest-isoform proteome: ${LONGEST}"
