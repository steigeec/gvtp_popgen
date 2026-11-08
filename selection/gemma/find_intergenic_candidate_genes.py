#!/usr/bin/env python3
"""
Collects every gene ID that ever appears as a flanking gene in a pure
snpEff "intergenic_region" annotation, across one or more outliers_*.ann.INFO
files, and resolves each to its GFF "gene" feature (scaffold, start, end).

intergenic_region names the nearest flanking gene at any distance, so these
genes are later kept only if the SNP is in LD with the gene body
(ld_scan_scaffold.py, curate_intergenic_by_ld.py).

The gene ID regex is the same as in Launch_find_outlier_annotate.sh.
Transcript-style IDs are resolved to gene IDs as in resolve_gene_ids.R;
tokens that cannot be resolved are dropped (logged to stderr).

Usage:
  python3 find_intergenic_candidate_genes.py outliers_*.ann.INFO \\
      --out candidate_gene_intervals.tsv
"""
import argparse
import re
import sys

BASE_DIR = "/private/groups/meyerlab/esteiger"

GFF = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"

GENE_TOKEN_RE = re.compile(
    r"MSTRG\.[0-9]+(?:\.p?[0-9]+)?|IDmodified-gene-[0-9]+|g[0-9]+|DEINVI_[0-9]+\.[0-9]+"
)


def normalize_mstrg(token):
    # matches the sed "s/\(MSTRG\.[0-9]\+\)\.[0-9]\+.*/\1/" truncation
    return re.sub(r"^(MSTRG\.[0-9]+)\.[0-9]+.*", r"\1", token)


def load_gff():
    """Returns (gene_intervals: {gene_id: (scaffold, start, end)},
                tx2gene: {transcript_id: gene_id})"""
    gene_intervals = {}
    tx2gene = {}
    with open(GFF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            scaffold, _, feature, start, end = fields[0], fields[1], fields[2], fields[3], fields[4]
            attrs = fields[8]
            if feature == "gene":
                m = re.search(r"ID=([^;]+)", attrs)
                if m:
                    gene_intervals[m.group(1)] = (scaffold, int(start), int(end))
            elif feature == "mRNA":
                pid = re.search(r"ID=([^;]+)", attrs)
                parent = re.search(r"Parent=([^;]+)", attrs)
                if pid and parent:
                    tx2gene[pid.group(1)] = parent.group(1)
    return gene_intervals, tx2gene


def resolve_to_gene(token, gene_intervals, tx2gene):
    """Same fallback chain as resolve_gene_ids.R: raw id, then common
    AUGUSTUS/transdecoder transcript suffixes, then a transcript-id prefix
    match for truncated ids like "DEINVI_00006.26586"."""
    token = normalize_mstrg(token)
    if token in gene_intervals:
        return token
    for suffix in (".t1", ".1"):
        cand = token + suffix
        if cand in tx2gene:
            return tx2gene[cand]
    prefix = token + "."
    for tx_id, gene_id in tx2gene.items():
        if tx_id.startswith(prefix):
            return gene_id
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ann_info_files", nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # (scaffold, raw_token) pairs pulled from intergenic_region sub-annotations
    raw_candidates = set()

    for path in args.ann_info_files:
        with open(path) as f:
            header = f.readline()
            if not header.startswith("CHROM"):
                f.seek(0)
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 5:
                    continue
                chrom, ann = cols[0], cols[4]
                for sub in ann.split(","):
                    parts = sub.split("|")
                    if len(parts) < 5 or parts[1] != "intergenic_region":
                        continue
                    for tok in GENE_TOKEN_RE.findall(parts[4]):
                        raw_candidates.add((chrom, tok))

    print(f"find_intergenic_candidate_genes.py: {len(raw_candidates)} raw "
          f"(scaffold, token) candidates from {len(args.ann_info_files)} file(s)",
          file=sys.stderr)

    gene_intervals, tx2gene = load_gff()

    resolved = {}  # gene_id -> (scaffold, start, end)
    unresolved = []
    for chrom, tok in sorted(raw_candidates):
        gene_id = resolve_to_gene(tok, gene_intervals, tx2gene)
        if gene_id is None:
            unresolved.append(tok)
            continue
        scaffold, start, end = gene_intervals[gene_id]
        if scaffold != chrom:
            # flanking gene must be on the same scaffold as the SNP
            unresolved.append(tok)
            continue
        resolved[gene_id] = (scaffold, start, end)

    if unresolved:
        preview = unresolved[:10]
        more = f" (+{len(unresolved) - 10} more)" if len(unresolved) > 10 else ""
        print(f"find_intergenic_candidate_genes.py: {len(unresolved)} token(s) "
              f"could not be resolved to a GFF gene, dropped: {preview}{more}",
              file=sys.stderr)

    with open(args.out, "w") as out:
        for gene_id, (scaffold, start, end) in sorted(resolved.items(), key=lambda kv: (kv[1][0], kv[1][1])):
            out.write(f"{scaffold}\t{gene_id}\t{start}\t{end}\n")

    print(f"find_intergenic_candidate_genes.py: wrote {len(resolved)} candidate "
          f"gene intervals to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
