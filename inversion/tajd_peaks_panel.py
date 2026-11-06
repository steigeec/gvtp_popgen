#!/usr/bin/env python3
"""
Stand-alone panel (candidate replacement for the bottom-right panel of
DEINVI_0006_stacked.py): the genome-wide distribution of Tajima's D in
Tranquillon Mtn (Strauss; 100-site windows from 4x_strauss.Diana.out, the
same track used by plot_TajD_extreme_combined.R), with each fog outlier peak
shown as a point below it at the median Tajima's D of the windows
overlapping the peak; the chr 6 fog peak (DEINVI_00006, GNS1/SNF1) is a
diamond.

Peaks: retained fog outlier SNPs (outliers_for_annot_40000_3_fog.txt)
grouped into peaks of SNPs <= PEAK_GAP apart (as in fog_deinvi6_figure.R).
A peak's windows = every window overlapping the span from its first to its
last outlier SNP (the nearest window if none overlaps, e.g. single-site
peaks).

Usage: tajd_peaks_panel.py
Output: OUT (figure), OUT_TSV (per-peak values)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAJD = "/private/groups/meyerlab/esteiger/selection/TajD/4x_strauss.Diana.out"
GEMMA = "/private/groups/meyerlab/esteiger/selection/gemma"
OUTLIERS = f"{GEMMA}/fog/outliers_for_annot_40000_3_fog.txt"
OUT = "/private/groups/meyerlab/esteiger/inversions/TajD_tranquillon_fog_peaks.jpg"
OUT_TSV = "/private/groups/meyerlab/esteiger/inversions/TajD_tranquillon_fog_peaks.tsv"
CHR6_PEAK = ("DEINVI_00006", 49_031_257, 49_146_391)
PEAK_GAP, BINS = 40_000, 120
LABEL_SIZE, TICK_SIZE = 24, 19                      # as in DEINVI_0006_stacked.py
rng = np.random.default_rng(1)


def main():
    td = pd.read_csv(TAJD, sep="\t", usecols=["chr", "w_start", "w_end", "TajimasD"]).dropna(subset=["TajimasD"])
    o = pd.read_csv(OUTLIERS, sep=" ", header=None).iloc[:, :3]; o.columns = ["chr", "snp", "pos"]
    smap = pd.read_csv(f"{GEMMA}/scaffold_snp_map.txt", sep=r"\s+", header=None, usecols=[0, 2]).drop_duplicates()
    o["scaffold"] = o.chr.map(dict(zip(smap[0], smap[2])))
    o = o.drop_duplicates("snp").sort_values(["scaffold", "pos"])
    o["peak"] = ((o.scaffold != o.scaffold.shift()) | (o.pos.diff() > PEAK_GAP)).cumsum()
    peaks = o.groupby("peak").agg(scaffold=("scaffold", "first"), start=("pos", "min"), end=("pos", "max"),
                                  n_snps=("pos", "size")).reset_index(drop=True)
    by_chr = {c: g for c, g in td.groupby("chr")}
    rows = []
    for p in peaks.itertuples():
        w = by_chr.get(p.scaffold)
        if w is None:
            continue
        hit = w[(w.w_end >= p.start) & (w.w_start <= p.end)]
        if not len(hit):
            mid = (w.w_start + w.w_end) / 2
            hit = w.iloc[[int(np.argmin(np.abs(mid - (p.start + p.end) / 2)))]]
        d = hit.TajimasD.median()
        rows.append(dict(scaffold=p.scaffold, start=p.start, end=p.end, n_snps=p.n_snps, n_windows=len(hit),
                         median_TajD=d, genome_percentile=(td.TajimasD < d).mean() * 100,
                         chr6_peak=(p.scaffold == CHR6_PEAK[0] and p.start == CHR6_PEAK[1])))
    r = pd.DataFrame(rows)
    r.to_csv(OUT_TSV, sep="\t", index=False)
    gm = td.TajimasD.median()
    print(f"{len(td):,} windows; genome median {gm:.3f}; {len(r)} peaks")
    print(r.sort_values("median_TajD").to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    lo, hi = np.percentile(td.TajimasD, [0.1, 99.9])
    ax.hist(td.TajimasD, bins=BINS, range=(lo, hi), density=True, color="#bdbdbd", edgecolor="none")   # outer 0.2% not drawn
    ax.axvline(gm, color="grey", ls="--", lw=1.5, zorder=1)
    top = ax.get_ylim()[1]
    other, c6 = r[~r.chr6_peak], r[r.chr6_peak]
    y0 = -0.10 * top                                       # strip of other peaks below the distribution
    ax.scatter(other.median_TajD, y0 + rng.uniform(-0.035, 0.035, len(other)) * top, s=110, marker="o",
               facecolor="white", edgecolor="black", linewidth=1.3, zorder=3)
    ax.scatter(c6.median_TajD, [y0 - 0.12 * top] * len(c6), s=260, marker="D", color="black", zorder=4)  # own row
    ax.set_ylim(y0 - 0.19 * top, top * 1.02)
    ax.set_yticks([t for t in ax.get_yticks() if 0 <= t <= top])
    ax.set_xlabel("Tajima's D", fontsize=LABEL_SIZE)
    ax.set_ylabel("Density", fontsize=LABEL_SIZE)
    ax.tick_params(labelsize=TICK_SIZE)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()
