#!/usr/bin/env python3
"""
Consistency of 3-group sample clustering along a scaffold.

Per winPCA window, samples are split into 3 groups by k-means on PC1 and PC2
(each scaled by sqrt(% variance explained)); groupings are compared between
windows with the adjusted Rand index (ARI).

Outputs (OUT_DIR):
  <tag>_cluster_consistency.jpg  window x window ARI heatmap, plus ARI between
                                 each window and the window 1 Mb further on
  <tag>_cluster_consistency.tsv  that 1-Mb-apart ARI profile

Usage: pc_cluster_consistency.py <winpca_prefix> <tag> [<tag2 prefix> <tag2> ...]
"""

import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

BASE_DIR = "/private/groups/meyerlab/esteiger"

OUT_DIR   = f"{BASE_DIR}/inversions/all_localities"
K         = 3
THIN      = 10      # use every 10th window (10 kb step, so 100 kb spacing) for the heatmap
LAG       = 100     # windows between compared pairs in the profile (100 x 10 kb = 1 Mb, non-overlapping)
# grey boxes per scaffold, as in DEINVI_0006_stacked.py (none drawn for others)
SHADE_REGIONS = {"DEINVI_00006": [(13.0, 29.4), (30.4, 43.7), (44.6, 49.3)]}


def window_groups(prefix):
    pc1 = pd.read_csv(f"{prefix}.pc_1.tsv.gz", sep="\t", index_col=0)
    pc2 = pd.read_csv(f"{prefix}.pc_2.tsv.gz", sep="\t", index_col=0)
    st = pd.read_csv(f"{prefix}.stat.tsv.gz", sep="\t", index_col=0)
    ok = pc1.notna().all(axis=1) & pc2.notna().all(axis=1)
    pc1, pc2, st = pc1[ok], pc2[ok], st[ok]
    w1 = np.sqrt(st["pc_1_ve"].values)[:, None]
    w2 = np.sqrt(st["pc_2_ve"].values)[:, None]
    X1, X2 = pc1.values * w1, pc2.values * w2
    labels = np.empty(X1.shape, dtype=int)
    for i in range(len(X1)):
        feats = np.column_stack([X1[i], X2[i]])
        labels[i] = KMeans(n_clusters=K, n_init=10, random_state=0).fit_predict(feats)
    return pc1.index.values, labels


def run(prefix, tag):
    print(f"[{tag}] clustering windows of {prefix} ...", flush=True)
    scaf = re.search(r"DEINVI_\d+", prefix).group(0)
    shade = SHADE_REGIONS.get(scaf, [])
    pos, labels = window_groups(prefix)
    mb = pos / 1e6

    # 1-Mb-apart profile
    prof = np.array([adjusted_rand_score(labels[i], labels[i + LAG]) for i in range(len(labels) - LAG)])
    prof_mid = (mb[:-LAG] + mb[LAG:]) / 2
    pd.DataFrame({"midpoint_Mb": prof_mid, "ARI_1Mb_apart": prof}).to_csv(
        f"{OUT_DIR}/{tag}_cluster_consistency.tsv", sep="\t", index=False)

    # window x window heatmap on thinned windows
    sub = labels[::THIN]
    sub_mb = mb[::THIN]
    n = len(sub)
    ari = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            ari[i, j] = ari[j, i] = adjusted_rand_score(sub[i], sub[j])

    fig, (ax_h, ax_p) = plt.subplots(2, 1, figsize=(10, 12.5), sharex=True,
                                     gridspec_kw={"height_ratios": [4, 1], "hspace": 0.05})
    ext = [sub_mb[0], sub_mb[-1], sub_mb[-1], sub_mb[0]]
    im = ax_h.imshow(ari, extent=ext, cmap="magma", vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    for lo, hi in shade:
        for v in (lo, hi):
            ax_h.axvline(v, color="white", lw=0.6, ls="--")
            ax_h.axhline(v, color="white", lw=0.6, ls="--")
    ax_h.set_ylabel("Position (Mb)", fontsize=14)
    ax_h.set_title(f"{tag}: agreement of 3-group splits (PC1+PC2) between windows (ARI)", fontsize=13)
    cb = fig.colorbar(im, ax=[ax_h, ax_p], fraction=0.03, pad=0.02)
    cb.set_label("Adjusted Rand index", fontsize=12)

    for lo, hi in shade:
        ax_p.axvspan(lo, hi, color="#808080", alpha=0.25, zorder=0)
    ax_p.plot(prof_mid, prof, color="black", lw=1)
    ax_p.set_ylim(-0.05, 1.05)
    ax_p.set_ylabel("ARI vs window\n1 Mb on", fontsize=12)
    ax_p.set_xlabel(f"Position on chromosome {int(scaf.split('_')[1])} (Mb)", fontsize=14)
    ax_p.set_xlim(sub_mb[0], sub_mb[-1])

    out = f"{OUT_DIR}/{tag}_cluster_consistency.jpg"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[{tag}] saved {out}")


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) < 2 or len(args) % 2:
        sys.exit(__doc__)
    for prefix, tag in zip(args[::2], args[1::2]):
        run(prefix, tag)
