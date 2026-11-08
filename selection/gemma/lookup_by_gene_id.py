#!/usr/bin/env python3
"""
Looks up per-gene hits (BLAST orthologs, GO annotation) by exact gene ID,
in place of `grep -Ff gene_list.txt file.tsv`.

By default each gene is mapped to its transcript ID via the GFF mRNA lines
(Parent=<gene>, ID=<transcript>), for reference files keyed by transcript ID
(e.g. BLAST hits). Gene IDs not found in the GFF are used as-is, with a warning.
Gene and transcript IDs often differ (e.g. gene g3013, transcript IDmodified-mRNA-287389).

--by-gene matches gene IDs directly, for files keyed by gene ID
(e.g. the combined GO annotation table).

Usage: python3 lookup_by_gene_id.py [--by-gene] <gene_list.txt> <reference.tsv> [key_column]
  key_column: 1-based column in reference.tsv holding the lookup key
              (default: 1)
Prints matching lines from reference.tsv to stdout.
"""
import sys
import re

BASE_DIR = "/private/groups/meyerlab/esteiger"

GFF = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"

def main():
    args = sys.argv[1:]
    by_gene = "--by-gene" in args
    args = [a for a in args if a != "--by-gene"]

    if len(args) < 2:
        sys.exit("Usage: lookup_by_gene_id.py [--by-gene] <gene_list.txt> <reference.tsv> [key_column]")
    gene_list_file, ref_file = args[0], args[1]
    key_col = int(args[2]) - 1 if len(args) > 2 else 0

    gene_ids = set()
    with open(gene_list_file) as f:
        for line in f:
            g = line.strip().rstrip(".")
            if g:
                gene_ids.add(g)

    if by_gene:
        target_transcripts = set(gene_ids)
    else:
        # gene_id -> transcript id from mRNA lines (the gene line's mRNA_id= attribute can be stale)
        mrna_real = {}
        with open(GFF) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "mRNA":
                    continue
                parent = re.search(r"Parent=([^;]+)", fields[8])
                mid = re.search(r"ID=([^;]+)", fields[8])
                if parent and mid:
                    mrna_real[parent.group(1)] = mid.group(1)

        target_transcripts = set()
        unresolved = []
        for g in gene_ids:
            t = mrna_real.get(g)
            if t:
                target_transcripts.add(t)
            else:
                unresolved.append(g)
                target_transcripts.add(g)  # fall back to the raw gene id

        if unresolved:
            preview = unresolved[:10]
            more = f" (+{len(unresolved) - 10} more)" if len(unresolved) > 10 else ""
            print(f"lookup_by_gene_id.py: {len(unresolved)} of {len(gene_ids)} gene(s) not found in GFF, "
                  f"falling back to raw id: {preview}{more}", file=sys.stderr)

    with open(ref_file) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= key_col:
                continue
            if cols[key_col] in target_transcripts:
                sys.stdout.write(line)

if __name__ == "__main__":
    main()
