#!/usr/bin/env python3
"""
Correct the mRNA_id= attribute on "gene" lines of
remove_duplicates_rm_internalstops_longest_isoform.gff, which for many genes
does not match the ID of the gene's own child mRNA line. Each gene has
exactly one mRNA; all other lines are copied unchanged.

Usage: python3 fix_gene_mRNA_id.py <input.gff> <output.gff>
"""
import sys
import re

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: fix_gene_mRNA_id.py <input.gff> <output.gff>")
    in_gff, out_gff = sys.argv[1], sys.argv[2]

    # Pass 1: collect the real mRNA ID for each gene, keyed by Parent=
    mrna_real = {}
    with open(in_gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "mRNA":
                continue
            attrs = fields[8]
            parent = re.search(r"Parent=([^;]+)", attrs)
            mid = re.search(r"ID=([^;]+)", attrs)
            if parent and mid:
                mrna_real[parent.group(1)] = mid.group(1)

    # Pass 2: rewrite gene lines' mRNA_id= to the real value; copy everything else verbatim
    n_fixed = 0
    n_genes = 0
    with open(in_gff) as f, open(out_gff, "w") as out:
        for line in f:
            fields = line.rstrip("\n").split("\t") if not line.startswith("#") else None
            if fields is None or len(fields) < 9 or fields[2] != "gene":
                out.write(line)
                continue

            n_genes += 1
            attrs = fields[8]
            gid = re.search(r"gene_id=([^;]+)", attrs)
            claimed = re.search(r"mRNA_id=([^;]+)", attrs)
            real = mrna_real.get(gid.group(1)) if gid else None

            if gid and claimed and real and claimed.group(1) != real:
                n_fixed += 1
                new_attrs = attrs[:claimed.start()] + f"mRNA_id={real}" + attrs[claimed.end():]
                fields[8] = new_attrs
                out.write("\t".join(fields) + "\n")
            else:
                out.write(line)

    print(f"Genes processed: {n_genes}")
    print(f"Gene lines fixed (mRNA_id corrected): {n_fixed}")
    print(f"Wrote corrected GFF to {out_gff}")

if __name__ == "__main__":
    main()
