#!/usr/bin/env python3
"""Map candidate genes in gene_summary_all.tsv to transcript IDs in the
longest-isoform annotation GFF. The annotation has one mRNA per gene, so a
transcript's TPM is the gene's TPM.

Input:  gene_summary_all.tsv (from selection/gemma/build_gene_summary_all.py),
        annotation GFF
Output: RNA_DIR/tissue_expression/results/candidate_gene_to_transcript.tsv
Usage:  python3 build_gene_map.py
"""
import re
import csv

RNA_DIR = "/scratch1/emma/RNA-seq"
GENOME_DIR = "/redser5/projects/GVTP/genome/gene_annotations/easel"

GFF = f"{GENOME_DIR}/fixes/longest_isoform.renamed.rm_internal_stops.gff"
SUMMARY = f"{RNA_DIR}/gene_summary_all.tsv"
OUT = f"{RNA_DIR}/tissue_expression/results/candidate_gene_to_transcript.tsv"

attr_re = re.compile(r'(\w+)=([^;]+)')


def parse_attrs(field):
    return dict(attr_re.findall(field))


# transcript_id -> (parent_gene_id, mrna_gene_tag)
tx_info = {}
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9 or cols[2] != "mRNA":
            continue
        attrs = parse_attrs(cols[8])
        tx_id = attrs.get("ID")
        parent = attrs.get("Parent")
        gene_tag = attrs.get("gene_id")
        if tx_id:
            tx_info[tx_id] = (parent, gene_tag)

# gene id -> transcript ids, by Parent and by gene_id tag
by_parent = {}
by_tag = {}
for tx_id, (parent, tag) in tx_info.items():
    by_parent.setdefault(parent, []).append(tx_id)
    by_tag.setdefault(tag, []).append(tx_id)


def full_id(name):
    return name if name.startswith("DEINVI_") else "DEINVI-" + name


rows = []
with open(SUMMARY) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for r in reader:
        name = r["Gene_Name"]
        if not name or name == "NA":
            continue
        gid = full_id(name)
        tx_ids = by_parent.get(gid) or by_tag.get(gid) or []
        match_source = "parent" if by_parent.get(gid) else ("mrna_gene_tag" if by_tag.get(gid) else "NONE")
        rows.append({
            "Variable": r["Variable"],
            "Gene_Name": name,
            "full_gff_id": gid,
            "match_source": match_source,
            "transcript_ids": ",".join(tx_ids) if tx_ids else "NA",
            "Acronym": r.get("Acronym", ""),
            "Functional_Description": r.get("Functional_Description", ""),
        })

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

n_total = len(rows)
n_matched = sum(1 for r in rows if r["transcript_ids"] != "NA")
n_multi = sum(1 for r in rows if r["transcript_ids"].count(",") > 0)
print(f"Total gene_summary rows processed: {n_total}")
print(f"Matched to >=1 transcript: {n_matched}")
print(f"Rows with >1 transcript match: {n_multi}")
print(f"Wrote: {OUT}")
