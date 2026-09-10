#!/usr/bin/env python3
"""
Drop-in replacement for the
  cut -f5 <ann.INFO> | grep -oE "MSTRG\\...." | sed "..."
step in Launch_find_outlier_annotate.sh. Prints raw gene ID tokens to
stdout, one per match, same format as before -- still meant to be piped
into resolve_gene_ids.R.

The difference: for a SNP whose ONLY snpEff annotation is
intergenic_region (no bounded upstream_gene_variant/downstream_gene_variant,
no genic annotation), the flanking gene(s) named by intergenic_region are
only emitted if find_intergenic_candidate_genes.py + ld_scan_scaffold.py
found the SNP's own position in that gene's LD-linked position set (i.e.
some SNP inside the gene body is in measurable LD with this exact outlier
SNP). Genic/bounded-window annotations are untouched -- extracted the same
way as before, unconditionally.

Requires ld_curation/candidate_gene_intervals.tsv and
ld_curation/linked_positions/<scaffold>.linked.tsv to already exist (see
find_intergenic_candidate_genes.py and launch_ld_scan.sh). If a needed
scaffold has no linked_positions file yet, every intergenic-only call on
that scaffold is dropped (logged to stderr) rather than silently kept.

Usage:
  python3 curate_intergenic_by_ld.py <ann.INFO file> \\
      --ld-curation-dir /path/to/ld_curation [--r2-threshold 0.7]
"""
import argparse
import re
import sys
from pathlib import Path

GFF = "/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"

GENE_TOKEN_RE = re.compile(
    r"MSTRG\.[0-9]+(?:\.p?[0-9]+)?|IDmodified-gene-[0-9]+|g[0-9]+|DEINVI_[0-9]+\.[0-9]+"
)


def normalize_mstrg(token):
    return re.sub(r"^(MSTRG\.[0-9]+)\.[0-9]+.*", r"\1", token)


def load_gff_resolver():
    gene_ids = set()
    tx2gene = {}
    with open(GFF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature, attrs = fields[2], fields[8]
            if feature == "gene":
                m = re.search(r"ID=([^;]+)", attrs)
                if m:
                    gene_ids.add(m.group(1))
            elif feature == "mRNA":
                pid = re.search(r"ID=([^;]+)", attrs)
                parent = re.search(r"Parent=([^;]+)", attrs)
                if pid and parent:
                    tx2gene[pid.group(1)] = parent.group(1)
    return gene_ids, tx2gene


def resolve_to_gene(token, gene_ids, tx2gene):
    token = normalize_mstrg(token)
    if token in gene_ids:
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


def load_linked_positions(ld_curation_dir, scaffolds, r2_threshold):
    """(scaffold, gene_id) -> set of linked positions with r2 >= threshold"""
    linked_dir = Path(ld_curation_dir) / "linked_positions"
    result = {}
    missing = []
    for scaffold in scaffolds:
        path = linked_dir / f"{scaffold}.linked.tsv"
        if not path.exists():
            missing.append(scaffold)
            continue
        with open(path) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 3:
                    continue
                gene_id, pos, r2 = cols[0], int(cols[1]), float(cols[2])
                if r2 < r2_threshold:
                    continue
                result.setdefault((scaffold, gene_id), set()).add(pos)
    if missing:
        print(f"curate_intergenic_by_ld.py: no LD scan output for {len(missing)} "
              f"scaffold(s), intergenic-only calls there will be dropped: "
              f"{missing[:10]}{'...' if len(missing) > 10 else ''}", file=sys.stderr)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ann_info_file")
    ap.add_argument("--ld-curation-dir", required=True)
    ap.add_argument("--r2-threshold", type=float, default=0.7)
    args = ap.parse_args()

    rows = []
    with open(args.ann_info_file) as f:
        header = f.readline()
        if not header.startswith("CHROM"):
            f.seek(0)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            rows.append((cols[0], int(cols[1]), cols[4]))

    scaffolds = sorted({r[0] for r in rows})
    linked = load_linked_positions(args.ld_curation_dir, scaffolds, args.r2_threshold)

    gene_ids_gff, tx2gene = load_gff_resolver()

    n_kept_intergenic = 0
    n_dropped_intergenic = 0
    out = sys.stdout
    for chrom, pos, ann in rows:
        subs = [s.split("|") for s in ann.split(",")]
        subs = [s for s in subs if len(s) >= 5]
        has_non_intergenic = any(s[1] != "intergenic_region" for s in subs)

        for s in subs:
            if s[1] != "intergenic_region":
                # Only the Gene_ID subfield (s[4]) -- not the whole entry
                # ("|".join(s)) -- since Feature_ID (s[6]) can carry a raw
                # AUGUSTUS transcript stem like "g1478.t2" for a gene that
                # was since renamed to "IDmodified-gene-N". Scanning the
                # whole entry matches that stem as a spurious extra "g1478"
                # gene call, which silently collides with a real, unrelated
                # gene elsewhere in the genome that kept id "g1478".
                for tok in GENE_TOKEN_RE.findall(s[4]):
                    out.write(tok + "\n")
                continue

            if has_non_intergenic:
                # redundant with a bounded annotation already emitted above
                continue

            for tok in GENE_TOKEN_RE.findall(s[4]):
                gene_id = resolve_to_gene(tok, gene_ids_gff, tx2gene)
                if gene_id is None:
                    n_dropped_intergenic += 1
                    continue
                positions = linked.get((chrom, gene_id))
                if positions and pos in positions:
                    out.write(gene_id + "\n")
                    n_kept_intergenic += 1
                else:
                    n_dropped_intergenic += 1

    print(f"curate_intergenic_by_ld.py: {n_kept_intergenic} intergenic-only gene "
          f"call(s) kept (LD-verified), {n_dropped_intergenic} dropped (no LD "
          f"evidence)", file=sys.stderr)


if __name__ == "__main__":
    main()
