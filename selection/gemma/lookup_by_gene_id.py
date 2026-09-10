#!/usr/bin/env python3
"""
Robust replacement for `grep -Ff gene_list.txt some_file.tsv` when looking up
per-gene hits (BLAST orthologs, GO annotation) by gene ID.

Plain substring grep is unsafe here: a gene's own ID is often NOT the same
string as its real transcript ID in these reference files (BLAST hits and
the GO annotation table are both keyed by transcript ID). For example gene
"g3013"'s real transcript is "IDmodified-mRNA-287389" — grep -Ff would never
match that gene at all, or worse, could cross-match an unrelated entry that
happens to share a substring. See fix_gene_mRNA_id.py for the same
underlying issue in the GFF itself.

By default this resolves each gene in gene_list.txt to its real canonical
transcript ID via the GFF's own "mRNA" feature lines (Parent=<gene> /
ID=<real transcript id>) before looking anything up, so it can't silently
miss or cross-match on a stale/wrong ID. Genes not found in the GFF fall
back to their raw ID (with a warning), matching this project's established
fallback convention. Use this mode for reference files keyed by transcript
ID, e.g. BLAST/orthology hits (run against CDS/transcript sequences).

Pass --by-gene to skip transcript resolution and match gene_list.txt
entries directly against the reference file's key column instead. Use this
mode for reference files already keyed by gene ID, e.g. the combined GO
annotation table (convert_eggnog_to_geneid_go.py / harvest_blast_go.py both
emit gene-id-keyed output) — resolving to a transcript ID first would mean
never matching anything in a gene-id-keyed file.

Usage: python3 lookup_by_gene_id.py [--by-gene] <gene_list.txt> <reference.tsv> [key_column]
  key_column: 1-based column in reference.tsv holding the lookup key
              (default: 1)
Prints matching lines from reference.tsv to stdout.
"""
import sys
import re

GFF = "/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"

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
        # gene_id -> real transcript id, from the GFF's own mRNA feature lines
        # (NOT the gene line's own mRNA_id= attribute, which can be stale)
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
