#!/usr/bin/env python3
"""
For one scaffold, streams its ngsLD pairwise output (DEINVI_XXXXX_prune.ld.gz,
columns: site1 site2 dist r2_ExpG D Dp r2 -- all pairs within 300kb, see
ngsLD.sh) exactly once and records, for each (variable, outlier SNP position)
seed on that scaffold, every position with r2 >= --r2-threshold to that seed.

Modeled directly on ld_scan_scaffold.py, but keyed by an exact seed position
(dict lookup) instead of gene-body interval containment (bisect). Seeds come
from outlier_snp_seeds.tsv (var_name, scaffold, pos), built by joining each
environmental variable's unfiltered_outliers_for_annot_<var>.txt against
scaffold_snp_map.txt. Kept separate per variable throughout (not pooled),
matching the rest of the Tajima's D pipeline's per-variable outlier
restriction.

Usage:
  python3 ld_scan_outlier_snps.py --scaffold DEINVI_00006 \\
      --ld-gz /path/DEINVI_00006_prune.ld.gz \\
      --seeds outlier_snp_seeds.tsv \\
      --r2-threshold 0.4 --threads 8 \\
      --out linked_snp_positions/DEINVI_00006.linked.tsv
"""
import argparse
import subprocess
import sys
from collections import defaultdict


def load_seeds(path, scaffold):
    """Returns {pos: set(var_names)} for seeds on this scaffold."""
    seeds = defaultdict(set)
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3 or cols[1] != scaffold:
                continue
            seeds[int(cols[2])].add(cols[0])
    return seeds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scaffold", required=True)
    ap.add_argument("--ld-gz", required=True)
    ap.add_argument("--seeds", required=True)
    ap.add_argument("--r2-threshold", type=float, default=0.4)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seeds = load_seeds(args.seeds, args.scaffold)
    if not seeds:
        print(f"ld_scan_outlier_snps.py: no outlier-SNP seeds on {args.scaffold}, "
              f"writing empty output", file=sys.stderr)
        open(args.out, "w").close()
        return

    # (var_name, seed_pos) -> {linked_pos: max_r2}
    linked = defaultdict(dict)

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

        vs1 = seeds.get(pos1)
        if vs1:
            for v in vs1:
                d = linked[(v, pos1)]
                if r2 > d.get(pos2, -1.0):
                    d[pos2] = r2

        vs2 = seeds.get(pos2)
        if vs2:
            for v in vs2:
                d = linked[(v, pos2)]
                if r2 > d.get(pos1, -1.0):
                    d[pos1] = r2

    ret = proc.wait()
    if ret != 0:
        sys.exit(f"ld_scan_outlier_snps.py: unpigz exited {ret} for {args.ld_gz}")

    with open(args.out, "w") as out:
        for (var_name, seed_pos) in sorted(linked):
            for pos, r2 in sorted(linked[(var_name, seed_pos)].items()):
                out.write(f"{var_name}\t{args.scaffold}\t{seed_pos}\t{pos}\t{r2:.4f}\n")

    n_linked = sum(len(d) for d in linked.values())
    print(f"ld_scan_outlier_snps.py: {args.scaffold} — {n_lines} LD pairs scanned, "
          f"{n_linked} linked positions across {len(seeds)} outlier-SNP seeds",
          file=sys.stderr)


if __name__ == "__main__":
    main()
