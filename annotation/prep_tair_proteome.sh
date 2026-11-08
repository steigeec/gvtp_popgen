#!/bin/bash
# Download the Arabidopsis TAIR10 proteome (Ensembl Plants release 63) and
# keep the longest protein per AGI gene, for OrthoFinder. Skips the download
# if the raw file already exists.
#
# Usage: bash prep_tair_proteome.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

OUTDIR=${BASE_DIR}/genome/annotation
RELEASE=63
URL="https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-${RELEASE}/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"

RAW_GZ="${OUTDIR}/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"
RAW="${OUTDIR}/Arabidopsis_thaliana.TAIR10.pep.all.fa"
LONGEST="${OUTDIR}/Arabidopsis_thaliana.TAIR10.pep.longest.fa"

if [ -s "$RAW_GZ" ]; then
    echo "Raw TAIR10 proteome already downloaded at ${RAW_GZ} — skipping download."
else
    echo "Downloading ${URL}..."
    curl -sS --fail -o "$RAW_GZ" "$URL"
fi

gunzip -k -f "$RAW_GZ"

# Headers: ">AT5G16970.1 pep ... gene:AT5G16970 ..."; keep the longest protein per gene
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

print(f"Wrote {len(best)} longest-isoform AT gene proteins to {out_file}")
PY

echo "Done. Longest-isoform Arabidopsis proteome: ${LONGEST}"
