#!/usr/bin/env python3
"""
Pick individuals with an unambiguous PC1/PC2 group in each candidate block on
DEINVI_00006, for the read-level breakpoint check (breakpoint_reads.sh /
breakpoint_summary.py).

For each block: samples are split into 3 groups (k-means on PC1 and PC2, each
scaled by sqrt(% variance explained), as in pc_cluster_consistency.py) in
every window of the block. The core window's split defines the groups; in
every other window, group labels are matched to the core's by maximum overlap
(Hungarian assignment), so PC1/PC2 swaps and label order don't matter. A
sample is "consistent" if it sits in its core group in >= MIN_CONSISTENCY of
the block's windows. Groups are ordered by mean PC1 at the core window: A
(low), Het (middle), B (high) -- which homozygous group matches the
reference assembly is unknown until the read check.

Per block and group, every consistent sample (or the TOP_N with the highest coverage
if TOP_N is set; Cov in env_var.txt) with an existing BAM is kept.

Output: OUT_DIR/breakpoint_samples.tsv (block, group, sample, consistency, cov)
"""

import os

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans

PREFIX   = "/private/groups/meyerlab/esteiger/inversions/all_localities/DEINVI_00006_GL_maf5_flipped"
ENV      = "/private/groups/meyerlab/esteiger/env_var.txt"
BAM_DIR  = "/private/groups/meyerlab/esteiger/mapping/merged"
OUT_DIR  = "/private/groups/meyerlab/esteiger/inversions/breakpoints"
MIN_CONSISTENCY = 0.9
TOP_N    = None       # max samples per block and group (None = every consistent sample)
STEP     = 5            # use every 5th window (50 kb spacing)

# block name: (start Mb, end Mb, core window Mb) -- spans from the
# three-group block extents (all-locality and Strauss-only runs agree)
BLOCKS = {
    "INV1_45.5-49.1":   (45.5, 49.1, 47.5),
    "blk_16.0-21.4":    (16.0, 21.4, 20.3),
    "blk_24.3-27.4":    (24.3, 27.4, 25.5),
    "blk_27.5-29.0":    (27.5, 29.0, 28.5),
    "blk_31.6-36.6":    (31.6, 36.6, 34.5),
}


def split(p1, p2, st, i):
    X = np.column_stack([p1.iloc[i] * np.sqrt(st.pc_1_ve.iloc[i]),
                         p2.iloc[i] * np.sqrt(st.pc_2_ve.iloc[i])])
    return KMeans(3, n_init=10, random_state=0).fit_predict(X)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    p1 = pd.read_csv(f"{PREFIX}.pc_1.tsv.gz", sep="\t", index_col=0)
    p2 = pd.read_csv(f"{PREFIX}.pc_2.tsv.gz", sep="\t", index_col=0)
    st = pd.read_csv(f"{PREFIX}.stat.tsv.gz", sep="\t", index_col=0)
    ok = p1.notna().all(axis=1) & p2.notna().all(axis=1)
    p1, p2, st = p1[ok], p2[ok], st[ok]
    mb = p1.index.values / 1e6
    samples = p1.columns

    env = pd.read_csv(ENV, sep="\t", dtype={"SampleID": str})
    cov = pd.to_numeric(env.set_index("SampleID")["Cov"], errors="coerce")
    cov = cov.groupby(level=0).max()
    sid = samples.str.replace("_rmdup.bam", "", regex=False)
    has_bam = np.array([os.path.exists(f"{BAM_DIR}/{s}") for s in samples])

    rows = []
    for block, (lo, hi, core) in BLOCKS.items():
        idx = np.where((mb >= lo) & (mb <= hi))[0][::STEP]
        ci = np.argmin(abs(mb - core))
        core_lab = split(p1, p2, st, ci)
        # name core groups A / Het / B by mean PC1
        order = np.argsort([p1.iloc[ci].values[core_lab == k].mean() for k in range(3)])
        name = {order[0]: "A", order[1]: "Het", order[2]: "B"}
        agree = np.zeros(len(samples))
        for i in idx:
            lab = split(p1, p2, st, i)
            overlap = np.array([[np.sum((lab == a) & (core_lab == b)) for b in range(3)] for a in range(3)])
            r, c = linear_sum_assignment(-overlap)
            mapped = np.empty_like(lab)
            for a, b in zip(r, c):
                mapped[lab == a] = b
            agree += (mapped == core_lab)
        frac = agree / len(idx)
        df = pd.DataFrame({"sample": samples, "group": [name[k] for k in core_lab],
                           "consistency": frac, "cov": cov.reindex(sid).values, "has_bam": has_bam})
        print(f"{block}: core-window group sizes " +
              ", ".join(f"{g}={n}" for g, n in df.group.value_counts().sort_index().items()) +
              "; consistent (>=90%): " +
              ", ".join(f"{g}={n}" for g, n in df[df.consistency >= MIN_CONSISTENCY].group.value_counts().sort_index().items()))
        for g in ["A", "Het", "B"]:
            pick = (df[(df.group == g) & (df.consistency >= MIN_CONSISTENCY) & df.has_bam & df["cov"].notna()]
                    .sort_values("cov", ascending=False).head(TOP_N if TOP_N else None))
            print(f"    {g:3s}: {len(pick)} picked, coverage {pick['cov'].min():.1f}-{pick['cov'].max():.1f}x "
                  f"(total {pick['cov'].sum():.0f}x)")
            for _, r in pick.iterrows():
                rows.append({"block": block, "group": g, "sample": r["sample"],
                             "consistency": round(r.consistency, 3), "cov": r["cov"]})
    out = pd.DataFrame(rows)
    out.to_csv(f"{OUT_DIR}/breakpoint_samples.tsv", sep="\t", index=False)
    print(f"\n{out['sample'].nunique()} unique samples across blocks -> {OUT_DIR}/breakpoint_samples.tsv")


if __name__ == "__main__":
    main()
