#!/usr/bin/env python3
"""
For one scaffold, streams its ngsLD pairwise output (DEINVI_XXXXX_prune.ld.gz,
columns: site1 site2 dist r2_ExpG D Dp r2, all pairs within 300kb -- see
ngsLD.sh) exactly once and records, for each candidate gene on that scaffold,
every position with r2 >= --r2-threshold to at least one SNP inside that
gene's body (start/end from candidate_gene_intervals.tsv, produced by
find_intergenic_candidate_genes.py).

This is deliberately decoupled from any specific outlier SNP or phenotype:
it answers "what positions does this gene's LD block actually reach", once,
genome-wide over the candidate set. curate_intergenic_by_ld.py then just
checks set membership per outlier SNP at pipeline time -- no rescanning.

Gene bodies are assumed non-overlapping on a given scaffold (true for this
GFF: longest-isoform-deduplicated AUGUSTUS/transdecoder models), so a
bisect over sorted start positions is used for containment instead of a
per-line linear scan -- required for speed, since the largest scaffold
files run to hundreds of millions of lines.

Usage:
  python3 ld_scan_scaffold.py --scaffold DEINVI_00021 \\
      --ld-gz /path/DEINVI_00021_prune.ld.gz \\
      --candidates candidate_gene_intervals.tsv \\
      --r2-threshold 0.7 --threads 8 \\
      --out linked_positions/DEINVI_00021.linked.tsv
"""
import argparse
import bisect
import subprocess
import sys


def load_candidates(path, scaffold):
    starts, ends, gene_ids = [], [], []
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4 or cols[0] != scaffold:
                continue
            starts.append(int(cols[2]))
            ends.append(int(cols[3]))
            gene_ids.append(cols[1])
    order = sorted(range(len(starts)), key=lambda i: starts[i])
    starts = [starts[i] for i in order]
    ends = [ends[i] for i in order]
    gene_ids = [gene_ids[i] for i in order]
    return starts, ends, gene_ids


def find_gene(pos, starts, ends, gene_ids):
    i = bisect.bisect_right(starts, pos) - 1
    if i >= 0 and pos <= ends[i]:
        return gene_ids[i]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scaffold", required=True)
    ap.add_argument("--ld-gz", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--r2-threshold", type=float, default=0.7)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    starts, ends, gene_ids = load_candidates(args.candidates, args.scaffold)
    if not starts:
        print(f"ld_scan_scaffold.py: no candidate genes on {args.scaffold}, "
              f"writing empty output", file=sys.stderr)
        open(args.out, "w").close()
        return

    # gene_id -> {position: max_r2}
    linked = {g: {} for g in set(gene_ids)}

    proc = subprocess.Popen(
        ["unpigz", "-p", str(args.threads), "-c", args.ld_gz],
        stdout=subprocess.PIPE, text=True, bufsize=1 << 20,
    )

    n_lines = 0
    first = True
    for line in proc.stdout:
        if first:
            first = False
            if line.startswith("site1"):
                continue
        n_lines += 1
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 7:
            continue
        try:
            pos1 = int(cols[0].rsplit(":", 1)[1])
            pos2 = int(cols[1].rsplit(":", 1)[1])
            r2 = float(cols[6])
        except (ValueError, IndexError):
            continue
        if r2 < args.r2_threshold:
            continue

        g1 = find_gene(pos1, starts, ends, gene_ids)
        if g1 is not None:
            d = linked[g1]
            if r2 > d.get(pos2, -1.0):
                d[pos2] = r2

        g2 = find_gene(pos2, starts, ends, gene_ids)
        if g2 is not None:
            d = linked[g2]
            if r2 > d.get(pos1, -1.0):
                d[pos1] = r2

    ret = proc.wait()
    if ret != 0:
        sys.exit(f"ld_scan_scaffold.py: unpigz exited {ret} for {args.ld_gz}")

    with open(args.out, "w") as out:
        for gene_id in sorted(linked):
            for pos, r2 in sorted(linked[gene_id].items()):
                out.write(f"{gene_id}\t{pos}\t{r2:.4f}\n")

    n_linked = sum(len(d) for d in linked.values())
    print(f"ld_scan_scaffold.py: {args.scaffold} — {n_lines} LD pairs scanned, "
          f"{n_linked} linked (gene,pos) pairs across {len(starts)} candidate genes",
          file=sys.stderr)


if __name__ == "__main__":
    main()
