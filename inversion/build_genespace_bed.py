#!/usr/bin/env python3
"""
Build a GENESPACE-format bed file (chr, start, end, id; 0-based start, no
header) from gff "gene" features, using a regex for the attribute that matches
the species' peptide fasta headers:
  - Ensembl Plants (Arabidopsis, Helianthus annuus, Lactuca, Solanum): 'ID=gene:([^;]+)'
  - NCBI Datasets (Helianthus anomalus): 'locus_tag=([^;]+)'
  - Deinandra: 'mRNA_id=([^;]+)' (peptide headers are transcript IDs)

Output is restricted to IDs present in the peptide fasta; mismatches in either
direction are reported, since GENESPACE would otherwise drop them silently.

Usage:
  python3 build_genespace_bed.py <gff_path> <id_regex> <pep_path> <out_bed>
"""
import sys
import re
import gzip

def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def main():
    if len(sys.argv) != 5:
        sys.exit("Usage: build_genespace_bed.py <gff_path> <id_regex> <pep_path> <out_bed>")
    gff_path, id_regex, pep_path, out_bed = sys.argv[1:5]
    id_pattern = re.compile(id_regex)

    pep_ids = set()
    with opener(pep_path) as f:
        for line in f:
            if line.startswith(">"):
                pep_ids.add(line[1:].split()[0].strip())

    rows = []
    seen_ids = set()
    with opener(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            m = id_pattern.search(fields[8])
            if not m:
                continue
            gene_id = m.group(1)
            seen_ids.add(gene_id)
            if gene_id not in pep_ids:
                continue
            chrom, start, end = fields[0], int(fields[3]) - 1, int(fields[4])
            rows.append((chrom, start, end, gene_id))

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_bed, "w") as out:
        for chrom, start, end, gene_id in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{gene_id}\n")

    missing_from_gff = pep_ids - seen_ids
    unmatched_gff = seen_ids - pep_ids
    print(f"Wrote {len(rows)} bed rows to {out_bed}")
    print(f"  peptide IDs with no gene-line match in gff at all: {len(missing_from_gff)}")
    print(f"  gff gene IDs with no matching peptide: {len(unmatched_gff)}")
    if missing_from_gff:
        print(f"  (sample missing-from-gff): {list(missing_from_gff)[:5]}")

if __name__ == "__main__":
    main()
