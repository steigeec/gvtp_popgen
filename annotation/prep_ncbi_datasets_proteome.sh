#!/bin/bash
# Download a proteome + GFF from NCBI Datasets (for species not on Ensembl
# Plants, e.g. Helianthus anomalus GCA_047456115.1) and keep the longest
# protein per gene, for OrthoFinder. Skips the download if the zip exists.
#
# Usage: bash prep_ncbi_datasets_proteome.sh <accession> <display_name>
#   accession:    NCBI Datasets genome accession, e.g. GCA_047456115.1
#   display_name: short name used for the output filename, e.g. "Helianthus_anomalus"
#
# Example:
#   bash prep_ncbi_datasets_proteome.sh GCA_047456115.1 Helianthus_anomalus

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

if [ $# -lt 2 ]; then
    echo "Usage: bash prep_ncbi_datasets_proteome.sh <accession> <display_name>" >&2
    exit 1
fi

ACCESSION="$1"
DISPLAY_NAME="$2"

OUTDIR=${BASE_DIR}/genome/annotation
WORKDIR="${OUTDIR}/ncbi_datasets_${DISPLAY_NAME}"
ZIP="${OUTDIR}/${DISPLAY_NAME}.ncbi_dataset.zip"
GFF_GZ="${OUTDIR}/${DISPLAY_NAME}.gff.gz"
RAW="${OUTDIR}/${DISPLAY_NAME}.pep.all.fa"
LONGEST="${OUTDIR}/${DISPLAY_NAME}.pep.longest.fa"

if [ -s "$ZIP" ]; then
    echo "NCBI Datasets zip for ${DISPLAY_NAME} already downloaded at ${ZIP} — skipping download."
else
    echo "Downloading ${ACCESSION} (protein FASTA + GFF) from NCBI Datasets..."
    curl -sS --fail --max-time 300 -L -o "$ZIP" \
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${ACCESSION}/download?include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_GFF"
fi

rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"
unzip -q -o "$ZIP" -d "$WORKDIR"

DATA_DIR=$(find "$WORKDIR/ncbi_dataset/data" -maxdepth 1 -type d -name "${ACCESSION}*")
PROTEIN_FAA="${DATA_DIR}/protein.faa"
GFF="${DATA_DIR}/genomic.gff"

if [ ! -s "$PROTEIN_FAA" ] || [ ! -s "$GFF" ]; then
    echo "ERROR: expected protein.faa and genomic.gff under ${DATA_DIR}" >&2
    exit 1
fi

cp "$PROTEIN_FAA" "$RAW"
gzip -kf "$GFF" -c > "$GFF_GZ"

# protein.faa headers carry no gene field; group isoforms by the GFF CDS locus_tag
python3 - "$GFF" "$RAW" "$LONGEST" <<'PY'
import sys
import re

gff_file, raw_file, out_file = sys.argv[1], sys.argv[2], sys.argv[3]

protein_to_gene = {}
with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        attrs = fields[8]
        pid_m = re.search(r"protein_id=([^;]+)", attrs)
        gene_m = re.search(r"locus_tag=([^;]+)", attrs)
        if pid_m and gene_m:
            protein_to_gene[pid_m.group(1)] = gene_m.group(1)

best = {}  # gene_id -> (length, protein_id, seq)
protein_id = None
seq_chunks = []

def flush():
    if protein_id is None:
        return
    gene_id = protein_to_gene.get(protein_id)
    if gene_id is None:
        return
    seq = "".join(seq_chunks)
    prev = best.get(gene_id)
    if prev is None or len(seq) > prev[0]:
        best[gene_id] = (len(seq), protein_id, seq)

with open(raw_file) as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            protein_id = line[1:].split()[0]
            seq_chunks = []
        else:
            seq_chunks.append(line)
    flush()

with open(out_file, "w") as out:
    for gene_id in sorted(best):
        _, protein_id, seq = best[gene_id]
        out.write(f">{gene_id}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i + 60] + "\n")

print(f"Wrote {len(best)} longest-isoform gene proteins to {out_file} "
      f"({len(protein_to_gene)} CDS->gene mappings parsed from GFF)")
PY

rm -rf "$WORKDIR"

echo "Done. Longest-isoform proteome: ${LONGEST}"
echo "GFF retained at: ${GFF_GZ}"
