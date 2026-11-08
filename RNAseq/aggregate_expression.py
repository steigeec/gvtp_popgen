#!/usr/bin/env python3
"""TPM of candidate genes from the per-sample salmon quant.sf files, as wide and
tidy tables with one row per gene_summary_all.tsv entry.

Input:  quant/<sample>/quant.sf (run_quant.sh),
        candidate_gene_to_transcript.tsv (build_gene_map.py), gene_summary_all.tsv
Output: RNA_DIR/tissue_expression/results/candidate_gene_TPM_{wide,tidy}.tsv
Usage:  python3 aggregate_expression.py
"""
import csv
import os

RNA_DIR = "/scratch1/emma/RNA-seq"
BASE = os.path.join(RNA_DIR, "tissue_expression")
QUANT_DIR = os.path.join(BASE, "quant")
MAP_FILE = os.path.join(BASE, "results/candidate_gene_to_transcript.tsv")
SUMMARY = os.path.join(RNA_DIR, "gene_summary_all.tsv")
OUT_WIDE = os.path.join(BASE, "results/candidate_gene_TPM_wide.tsv")
OUT_TIDY = os.path.join(BASE, "results/candidate_gene_TPM_tidy.tsv")

SAMPLE_META = {
    "GVTP-1Root_S2":  {"individual": "1", "tissue": "Root"},
    "GVTP-1Shoot_S1": {"individual": "1", "tissue": "Shoot"},
    "GVTP-2Root_S4":  {"individual": "2", "tissue": "Root"},
    "GVTP-2Shoot_S3": {"individual": "2", "tissue": "Shoot"},
    "GVTP-5Root_S6":  {"individual": "5", "tissue": "Root"},
    "GVTP-5Shoot_S5": {"individual": "5", "tissue": "Shoot"},
    "GVTPMLeaf_S7":   {"individual": "M", "tissue": "Leaf"},
    "GVTPMRoot_S8":   {"individual": "M", "tissue": "Root"},
}

# transcript_id -> {sample: tpm}
tx_tpm = {}
for sample in SAMPLE_META:
    sf = os.path.join(QUANT_DIR, sample, "quant.sf")
    with open(sf) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tx_tpm.setdefault(row["Name"], {})[sample] = float(row["TPM"])

# salmon's index keeps one copy of identical transcripts (e.g. conserved
# chloroplast ribosomal protein genes); map each dropped id to the kept one.
DUP_CLUSTERS = os.path.join(BASE, "ref/salmon_index/duplicate_clusters.tsv")
dup_to_retained = {}
with open(DUP_CLUSTERS) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for r in reader:
        dup_to_retained[r["DuplicateRef"]] = r["RetainedRef"]


def lookup_tpm(tx_id):
    if tx_id in tx_tpm:
        return tx_tpm[tx_id]
    if tx_id in dup_to_retained:
        return tx_tpm.get(dup_to_retained[tx_id])
    return None

# gene_summary rows (a gene can appear once per environmental variable)
summary_rows = []
with open(SUMMARY) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for r in reader:
        summary_rows.append(r)

gene_map = {}
with open(MAP_FILE) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for r in reader:
        gene_map[(r["Variable"], r["Gene_Name"])] = r["transcript_ids"]

samples = list(SAMPLE_META.keys())
wide_rows = []
tidy_rows = []
missing = []
for r in summary_rows:
    key = (r["Variable"], r["Gene_Name"])
    tx_ids_field = gene_map.get(key, "NA")
    if tx_ids_field in (None, "NA"):
        missing.append(key)
        continue
    tx_ids = tx_ids_field.split(",")
    # sum TPM across matched transcripts (normally just 1)
    sample_tpm = {s: 0.0 for s in samples}
    found_any = False
    for tx in tx_ids:
        d = lookup_tpm(tx)
        if d:
            found_any = True
            for s in samples:
                sample_tpm[s] += d.get(s, 0.0)
    if not found_any:
        missing.append(key)
        continue
    row = {
        "Variable": r["Variable"],
        "Gene_Name": r["Gene_Name"],
        "Acronym": r.get("Acronym", ""),
        "Functional_Description": r.get("Functional_Description", ""),
        "transcript_ids": tx_ids_field,
    }
    row.update({s: round(sample_tpm[s], 3) for s in samples})
    wide_rows.append(row)
    for s in samples:
        tidy_rows.append({
            "Variable": r["Variable"],
            "Gene_Name": r["Gene_Name"],
            "Acronym": r.get("Acronym", ""),
            "sample": s,
            "individual": SAMPLE_META[s]["individual"],
            "tissue": SAMPLE_META[s]["tissue"],
            "TPM": round(sample_tpm[s], 3),
        })

fieldnames = ["Variable", "Gene_Name", "Acronym", "Functional_Description", "transcript_ids"] + samples
with open(OUT_WIDE, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(wide_rows)

with open(OUT_TIDY, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["Variable", "Gene_Name", "Acronym", "sample", "individual", "tissue", "TPM"], delimiter="\t")
    w.writeheader()
    w.writerows(tidy_rows)

print(f"Wide table rows: {len(wide_rows)} -> {OUT_WIDE}")
print(f"Tidy table rows: {len(tidy_rows)} -> {OUT_TIDY}")
print(f"Missing/unquantified rows: {len(missing)}")
if missing:
    print(missing[:20])
