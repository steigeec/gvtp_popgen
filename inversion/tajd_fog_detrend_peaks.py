#!/usr/bin/env python3
"""
Stand-alone figure (formerly the bottom-right panel of DEINVI_0006_stacked.py):
Delta Tajima's D of the fog outlier peaks against the detrended SD of fog, per
locality, with one box per locality summarizing all fog peaks.

Delta Tajima's D of a peak in a population = median Tajima's D of the
100-site windows (4x_<pop>.Diana.out, as in plot_TajD_extreme_combined.R)
overlapping the peak's span (first to last outlier SNP; the nearest window if
none overlaps, e.g. single-SNP peaks) minus that population's genome-wide
median window Tajima's D. Peaks = retained fog outlier SNPs
(outliers_for_annot_40000_3_fog.txt) grouped at <= PEAK_GAP, as elsewhere.
x = mean detrended SD of fog over the population's Tajima's D sample (the
mean_fog_detrend column of plot_TajD_extreme_combined.R's FOG_X table).
OUT: box = median and quartiles of the peaks (whiskers 1.5 IQR); points = the
individual peaks.
OUT_NORM (reaction norms): each peak's five site values joined by a line, in
the style of DEINVI_0006_stacked.py panels G and I -- the chr 6 fog peak
(CHR6_PEAK) as diamonds and a black dashed line, the other peaks as points and
thin lines at OTHER_ALPHA; points coloured by site; black bars = site medians
across peaks.
OUT_TREND: the same reaction norms against the MAGNITUDE of long-term fog
change instead (prediction: sites whose fog has changed more, in either
direction, have more negative Delta D at fog outliers): |mean Kling fog trend
1950-2014| (gvtp_popgen/env/fog_trend.R; % points per decade) over the plants
in the population's Tajima's D sample (bamlist_4x_<pop>.txt), untrustworthy-
coordinate plants excluded.
OUT_FOG: the same reaction norms against mean long-term fog (mean_fog of the
FOG_X table: MODIS June-September 2000-2022 days/month, averaged over the
population's Tajima's D sample).
OUT_CWD: the same reaction norms against long-term climatic water deficit
(BCM v8 mean annual CWD, WY2000-2022, mm; env/BCM_CWD_extracted.txt,
extract_BCM_CWD.R), averaged over the same plants as fog change.

Usage: tajd_fog_detrend_peaks.py
Output: OUT, OUT_NORM, OUT_TREND, OUT_FOG, OUT_CWD (figures), OUT_TSV (per-peak values)
"""

import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, "/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma")
import temporal_null as tn

TAJD_DIR = "/private/groups/meyerlab/esteiger/selection/TajD"
GEMMA = "/private/groups/meyerlab/esteiger/selection/gemma"
OUTLIERS = f"{GEMMA}/fog/outliers_for_annot_40000_3_fog.txt"
FOG_X = f"{GEMMA}/TajD_delta_vs_fogDetrend_fogPeak.tsv"
OUT = "/private/groups/meyerlab/esteiger/inversions/TajD_fog_detrend_peaks.jpg"
OUT_NORM = "/private/groups/meyerlab/esteiger/inversions/TajD_fog_detrend_peaks_norms.jpg"
OUT_TREND = "/private/groups/meyerlab/esteiger/inversions/TajD_fog_change_peaks_norms.jpg"
OUT_FOG = "/private/groups/meyerlab/esteiger/inversions/TajD_fog_mean_peaks_norms.jpg"
OUT_CWD = "/private/groups/meyerlab/esteiger/inversions/TajD_cwd_peaks_norms.jpg"
CWD = "/private/groups/meyerlab/esteiger/env/BCM_CWD_extracted.txt"
FOG_TREND = "/private/groups/meyerlab/esteiger/env/fog_trend_by_sample.tsv"
BAMLIST_DIR = "/private/groups/meyerlab/esteiger/mapping/merged"
OUT_TSV = "/private/groups/meyerlab/esteiger/inversions/TajD_fog_detrend_peaks.tsv"
POPS = {"CDFW": "cdfw", "Strauss": "strauss", "Point_Conception": "conception",
        "VAFB_Lion_Head": "lionhead", "VAFB_Arguello": "arguello"}
PEAK_GAP, BOX_WIDTH = 40_000, 0.04
CHR6_PEAK = ("DEINVI_00006", 49_031_257)
OTHER_ALPHA, OTHER_LINE = 0.45, "#7f7f7f"           # as the non-chr 6 peaks in DEINVI_0006_stacked.py G/I
LABEL_SIZE, TICK_SIZE = 24, 19                      # as in DEINVI_0006_stacked.py
rng = np.random.default_rng(1)


def fog_peaks():
    o = pd.read_csv(OUTLIERS, sep=" ", header=None).iloc[:, :3]; o.columns = ["chr", "snp", "pos"]
    smap = pd.read_csv(f"{GEMMA}/scaffold_snp_map.txt", sep=r"\s+", header=None, usecols=[0, 2]).drop_duplicates()
    o["scaffold"] = o.chr.map(dict(zip(smap[0], smap[2])))
    o = o.drop_duplicates("snp").sort_values(["scaffold", "pos"])
    o["peak"] = ((o.scaffold != o.scaffold.shift()) | (o.pos.diff() > PEAK_GAP)).cumsum()
    return o.groupby("peak").agg(scaffold=("scaffold", "first"), start=("pos", "min"), end=("pos", "max")).reset_index(drop=True)


def main():
    peaks = fog_peaks()
    fx = pd.read_csv(FOG_X, sep="\t").drop_duplicates("population").set_index("population")
    trend = pd.read_csv(FOG_TREND, sep="\t", dtype={"SampleID": str}).set_index("SampleID")
    cwd = pd.read_csv(CWD, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    fog_change, cwd_mean = {}, {}
    for pop, code in POPS.items():
        ids = [l.strip().replace("_rmdup.bam", "") for l in open(f"{BAMLIST_DIR}/bamlist_4x_{code}.txt") if l.strip()]
        ids = [i for i in ids if i in trend.index and i not in tn.BAD_COORD_IDS]
        cwd_mean[pop] = cwd.loc[[i for i in ids if i in cwd.index], "BCM_CWD_mm_30y"].mean()
        fog_change[pop] = abs(trend.loc[ids, "kling_fog_trend_decade"].mean() * 100)   # |% points / decade|
    rows = []
    for pop, code in POPS.items():
        td = pd.read_csv(f"{TAJD_DIR}/4x_{code}.Diana.out", sep="\t", usecols=["chr", "w_start", "w_end", "TajimasD"]).dropna()
        bg = td.TajimasD.median()
        by_chr = {c: g for c, g in td.groupby("chr")}
        for p in peaks.itertuples():
            w = by_chr.get(p.scaffold)
            if w is None:
                continue
            hit = w[(w.w_end >= p.start) & (w.w_start <= p.end)]
            if not len(hit):
                mid = (w.w_start + w.w_end) / 2
                hit = w.iloc[[int(np.argmin(np.abs(mid - (p.start + p.end) / 2)))]]
            rows.append(dict(population=pop, site=fx.loc[pop, "site"], color=fx.loc[pop, "color"],
                             fog_detrend=fx.loc[pop, "mean_fog_detrend"], fog_change=fog_change[pop], fog_mean=fx.loc[pop, "mean_fog"], cwd=cwd_mean[pop], scaffold=p.scaffold, start=p.start,
                             end=p.end, n_windows=len(hit), delta_TajD=hit.TajimasD.median() - bg))
    r = pd.DataFrame(rows)
    r.to_csv(OUT_TSV, sep="\t", index=False)
    s = r.groupby("site").agg(fog_detrend=("fog_detrend", "first"), fog_change=("fog_change", "first"), fog_mean=("fog_mean", "first"), cwd=("cwd", "first"), n_peaks=("delta_TajD", "size"),
                              median=("delta_TajD", "median"), mean=("delta_TajD", "mean")).sort_values("fog_detrend")
    print(s.round(3).to_string())
    rho = pd.Series(s["median"].values).corr(pd.Series(s.fog_detrend.values), method="spearman")
    print(f"Spearman rho (site median vs detrended SD of fog, 5 sites) = {rho:.2f}")

    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.axhline(0, color="grey", linestyle="--", linewidth=1.5, zorder=0)
    for site, g in r.groupby("site"):
        x = g.fog_detrend.iloc[0]
        ax.boxplot([g.delta_TajD], positions=[x], widths=BOX_WIDTH, showfliers=False, patch_artist=True, manage_ticks=False,
                   boxprops=dict(facecolor=g.color.iloc[0], edgecolor="black", linewidth=1.5),
                   medianprops=dict(color="black", linewidth=2.5), whiskerprops=dict(linewidth=1.5),
                   capprops=dict(linewidth=0), zorder=2)
        ax.scatter(x + rng.uniform(-BOX_WIDTH / 3, BOX_WIDTH / 3, len(g)), g.delta_TajD, s=18, color="black",
                   alpha=0.5, edgecolor="none", zorder=3)
    xs = s.fog_detrend.values
    ax.set_xlim(xs.min() - 0.08, xs.max() + 0.08)
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    ax.set_xlabel("Detrended SD of fog (days/month)", fontsize=LABEL_SIZE)
    ax.set_ylabel("Δ Tajima's D", fontsize=LABEL_SIZE)
    ax.tick_params(labelsize=TICK_SIZE)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"Saved: {OUT}")

    # ---- reaction norms ----
    r["pk"] = list(zip(r.scaffold, r.start))
    norm_plot(r, s, "fog_detrend", "Detrended SD of fog (days/month)", 0.1, OUT_NORM)
    norm_plot(r, s, "fog_change", "Magnitude of fog change\n(% points/decade, 1950-2014)", None, OUT_TREND)
    norm_plot(r, s, "fog_mean", "Long-term fog (days/month)", None, OUT_FOG)
    norm_plot(r, s, "cwd", "Climatic water deficit (mm/year)", None, OUT_CWD)
    for col in ("fog_detrend", "fog_change", "fog_mean", "cwd"):
        sl = r.groupby("pk").apply(lambda g: np.polyfit(g[col], g.delta_TajD, 1)[0])
        rho = stats.spearmanr(s[col], s["median"])
        print(f"{col}: site medians Spearman rho = {rho.statistic:.2f} (p = {rho.pvalue:.2f}); per-peak slopes: "
              f"{(sl < 0).sum()} of {len(sl)} negative, median {sl.median():+.3f}; chr 6 slope {sl[CHR6_PEAK]:+.3f}; "
              f"Wilcoxon (slope < 0) p = {stats.wilcoxon(sl, alternative='less').pvalue:.3f}")


def norm_plot(r, s, col, xlabel, tick, out):
    """reaction norms: each peak's site values joined by a line; chr 6 = diamonds + dashed line"""
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.axhline(0, color="grey", linestyle="--", linewidth=1.5, zorder=0)
    for pk, g in r.sort_values(col).groupby("pk"):
        if pk == CHR6_PEAK:
            continue
        ax.plot(g[col], g.delta_TajD, color=OTHER_LINE, lw=0.9, alpha=OTHER_ALPHA, zorder=1)
        ax.scatter(g[col], g.delta_TajD, s=30, c=g.color, edgecolor="black", linewidth=0.5,
                   alpha=OTHER_ALPHA, zorder=2)
    c6 = r[r.pk == CHR6_PEAK].sort_values(col)
    ax.plot(c6[col], c6.delta_TajD, color="black", lw=2.5, ls="--", zorder=4)
    ax.scatter(c6[col], c6.delta_TajD, s=170, marker="D", c=c6.color, edgecolor="black",
               linewidth=1.3, zorder=5)
    ax.scatter(s[col], s["median"], s=900, marker="_", color="black", linewidth=3.5, zorder=3)
    xs = s[col].values
    pad = (xs.max() - xs.min()) * 0.12
    ax.set_xlim(xs.min() - pad, xs.max() + pad)
    if tick:
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(tick))
    ax.set_xlabel(xlabel, fontsize=LABEL_SIZE)
    ax.set_ylabel("\u0394 Tajima's D", fontsize=LABEL_SIZE)
    ax.tick_params(labelsize=TICK_SIZE)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
