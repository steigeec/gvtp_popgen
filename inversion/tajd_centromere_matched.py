#!/usr/bin/env python3
"""
Is the chr 6 fog peak's Tajima's D unusually low for its genomic setting?
(A test of the meiotic-drive interpretation: a driven haplotype sweeping
through populations should leave lower Tajima's D than comparable regions,
in every population.)

Peak statistic: per population, median Tajima's D of the 100-site windows
(4x_<pop>.Diana.out, as in plot_TajD_extreme_combined.R) overlapping the chr 6
peak (PEAK; 19 retained fog SNPs, GNS1/SNF1).
Centromere-matched null: N_BLOCKS random blocks of the same span placed
outside the centromere call but within the same distance band from the
centromere edge as the peak (temporal_null.cen_band of the peak's distance:
0-1.59 Mb), on scaffolds with a confident centromere call
(call_centromeres.py), excluding DEINVI_00006 and +-EXCLUDE_FLANK around
every fog outlier peak; blocks with < MIN_WINDOWS windows dropped. The same
block positions are used in every population.
Genome-wide null: the same number of random blocks anywhere on those
scaffolds (same exclusions), to show how much is explained by centromere
proximity alone.
Context: every other fog peak's median Tajima's D (windows overlapping its
SNP span), per population.
One-sided p (drive predicts lower): fraction of null blocks with median
Tajima's D <= the peak's.

Figure (OUT_FIG): each locus sits next to its appropriate null -- the other
fog peaks (mostly on chromosome arms) beside the genome-wide random blocks,
the chr 6 peak beside the centromere-matched blocks.
OUT_FIG_KDE: per population, a smoothed (Gaussian KDE) frequency distribution
of all genome-wide 100-site window Tajima's D values, with each fog peak's
median Tajima's D as a point below it (chr 6 = diamond). Note the scale
difference: a peak value is a median over many windows, so it is less
dispersed than single windows.

Usage: tajd_centromere_matched.py
Output: OUT_FIG, OUT_FIG_KDE, OUT_TSV
"""

import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

sys.path.insert(0, "/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma")
import temporal_null as tn

TAJD_DIR = "/private/groups/meyerlab/esteiger/selection/TajD"
GEMMA = "/private/groups/meyerlab/esteiger/selection/gemma"
OUTLIERS = f"{GEMMA}/fog/outliers_for_annot_40000_3_fog.txt"
OUT_FIG = "/private/groups/meyerlab/esteiger/inversions/centromere/TajD_centromere_matched.jpg"
OUT_FIG_KDE = "/private/groups/meyerlab/esteiger/inversions/centromere/TajD_genome_kde_fog_peaks.jpg"
OUT_TSV = "/private/groups/meyerlab/esteiger/inversions/centromere/TajD_centromere_matched.tsv"
POP_COLOR = {"Tranquillon Mtn": "#009e73", "Gaviota": "#e69f00", "Pt. Conception": "#0072b2",
             "Lion's Head": "#f0e442", "Pt. Arguello": "#d55e00"}
KDE_SAMPLE, OTHER_ALPHA = 50_000, 0.45
# plotting order, north-west to south-east (Dangermond has no Tajima's D track: 1 plant)
POPS = {"lionhead": "Lion's Head", "arguello": "Pt. Arguello", "strauss": "Tranquillon Mtn",
        "conception": "Pt. Conception", "cdfw": "Gaviota"}
PEAK_SIZE = 115               # chr 6 diamond marker area
PEAK = ("DEINVI_00006", 49_031_257, 49_146_391)
N_BLOCKS, MIN_WINDOWS, EXCLUDE_FLANK, PEAK_GAP, SEED = 3000, 20, 500_000, 40_000, 7
LABEL_SIZE, TICK_SIZE = 22, 17


def fog_peaks():
    o = pd.read_csv(OUTLIERS, sep=" ", header=None).iloc[:, :3]; o.columns = ["chr", "snp", "pos"]
    smap = pd.read_csv(f"{GEMMA}/scaffold_snp_map.txt", sep=r"\s+", header=None, usecols=[0, 2]).drop_duplicates()
    o["scaffold"] = o.chr.map(dict(zip(smap[0], smap[2])))
    o = o.drop_duplicates("snp").sort_values(["scaffold", "pos"])
    o["peak"] = ((o.scaffold != o.scaffold.shift()) | (o.pos.diff() > PEAK_GAP)).cumsum()
    return o.groupby("peak").agg(scaffold=("scaffold", "first"), start=("pos", "min"), end=("pos", "max")).reset_index(drop=True)


def draw_blocks(cen, span, peaks, lo, hi, matched, rng):
    """random block starts: within [lo, hi] of the centromere edge (outside the call) if matched, else anywhere"""
    iv = []
    for chrom, c in cen.iterrows():
        if chrom == PEAK[0]:
            continue
        if matched:
            cand = [(c.cen_start - hi - span, c.cen_start - lo - span), (c.cen_end + lo, c.cen_end + hi)]
        else:
            cand = [(0, c.len - span)]
        iv += [(chrom, max(0, a), min(c.len - span, b)) for a, b in cand if min(c.len - span, b) - max(0, a) > 0]
    w = np.array([b - a for _, a, b in iv], float)
    out = []
    while len(out) < N_BLOCKS * 3 and len(out) < 50 * N_BLOCKS:
        ch, a, b = iv[rng.choice(len(iv), p=w / w.sum())]
        s = int(rng.integers(a, b))
        near = peaks[(peaks.scaffold == ch) & (peaks.end + EXCLUDE_FLANK > s) & (peaks.start - EXCLUDE_FLANK < s + span)]
        if not len(near):
            out.append((ch, s, s + span))
    return pd.DataFrame(out, columns=["chr", "start", "end"])


def block_medians(td_by_chr, blocks):
    med, nwin = [], []
    for ch, s, e in blocks.itertuples(index=False):
        w = td_by_chr.get(ch)
        if w is None:
            med.append(np.nan); nwin.append(0); continue
        i0, i1 = np.searchsorted(w["w_end"].values, s), np.searchsorted(w["w_start"].values, e, side="right")
        v = w["TajimasD"].values[i0:i1]
        med.append(np.median(v) if len(v) else np.nan); nwin.append(len(v))
    return np.array(med), np.array(nwin)


def main():
    rng = np.random.default_rng(SEED)
    cen = tn.centromeres()
    peaks = fog_peaks()
    span = PEAK[2] - PEAK[1]
    dist = tn.cen_dist(PEAK[0], (PEAK[1] + PEAK[2]) / 2, cen)
    lo, hi = tn.cen_band(dist)
    lo = max(lo, 1)                                   # outside the centromere call, like the peak
    matched = draw_blocks(cen, span, peaks, lo, hi, True, rng)
    genome = draw_blocks(cen, span, peaks, 0, 0, False, rng)
    print(f"peak {dist / 1e6:.2f} Mb from its centromere; matched band {lo / 1e6:.2f}-{hi / 1e6:.2f} Mb; "
          f"block span {span / 1e3:.0f} kb; candidate blocks: {len(matched)} matched, {len(genome)} genome-wide")

    rows, plot = [], {}
    for code, name in POPS.items():
        td = pd.read_csv(f"{TAJD_DIR}/4x_{code}.Diana.out", sep="\t", usecols=["chr", "w_start", "w_end", "TajimasD"]).dropna()
        td = td.sort_values(["chr", "w_start"])
        by_chr = {c: g.reset_index(drop=True) for c, g in td.groupby("chr")}
        pm, pn = block_medians(by_chr, pd.DataFrame([PEAK], columns=["chr", "start", "end"]))
        mm, mn = block_medians(by_chr, matched); gm, gn = block_medians(by_chr, genome)
        mm, gm = mm[mn >= MIN_WINDOWS][:N_BLOCKS], gm[gn >= MIN_WINDOWS][:N_BLOCKS]
        other = []
        for p in peaks.itertuples():
            if (p.scaffold, p.start) == (PEAK[0], PEAK[1]):
                continue
            s, e = (p.start, p.end) if p.end - p.start >= 1000 else (p.start - 500, p.end + 500)
            m, n = block_medians(by_chr, pd.DataFrame([(p.scaffold, s, e)], columns=["chr", "start", "end"]))
            if n[0]:
                other.append(m[0])
        peak = pm[0]
        rows.append(dict(population=name, peak_windows=int(pn[0]), peak_median_D=peak, genome_window_median=td.TajimasD.median(),
                         n_matched_blocks=len(mm), matched_block_median=np.median(mm),
                         peak_percentile_matched=(mm < peak).mean() * 100, p_matched=(np.sum(mm <= peak) + 1) / (len(mm) + 1),
                         n_genome_blocks=len(gm), genome_block_median=np.median(gm),
                         peak_percentile_genome=(gm < peak).mean() * 100, p_genome=(np.sum(gm <= peak) + 1) / (len(gm) + 1),
                         other_fog_peaks_median=np.median(other),
                         peak_rank_among_fog_peaks=f"{int(np.sum(np.array(other) < peak)) + 1} of {len(other) + 1}"))
        plot[name] = dict(matched=mm, genome=gm, other=np.array(other), peak=peak,
                          windows=rng.choice(td.TajimasD.values, min(KDE_SAMPLE, len(td)), replace=False))
        print(rows[-1], flush=True)
    r = pd.DataFrame(rows)
    r.to_csv(OUT_TSV, sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(r.to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    fig, ax = plt.subplots(figsize=(12, 6))
    names = list(plot)
    x = np.arange(len(names))
    bp1 = ax.boxplot([plot[n]["genome"] for n in names], positions=x - 0.2, widths=0.3, showfliers=False,
                     patch_artist=True, boxprops=dict(facecolor="white", edgecolor="grey"),
                     medianprops=dict(color="grey"), whiskerprops=dict(color="grey"), capprops=dict(color="grey"))
    bp2 = ax.boxplot([plot[n]["matched"] for n in names], positions=x + 0.2, widths=0.3, showfliers=False,
                     patch_artist=True, boxprops=dict(facecolor="#bdbdbd", edgecolor="black"),
                     medianprops=dict(color="black", linewidth=2))
    for i, n in enumerate(names):
        o = plot[n]["other"]
        ax.scatter(i - 0.2 + rng.uniform(-0.08, 0.08, len(o)), o, s=28, facecolor=POP_COLOR[n],
                   edgecolor="black", linewidth=0.6, zorder=3)
        ax.scatter(i + 0.2, plot[n]["peak"], s=PEAK_SIZE, marker="D", facecolor=POP_COLOR[n],
                   edgecolor="black", linewidth=1.3, zorder=4)
    # legend proxies in neutral grey (fills follow each population's colour)
    h_other = plt.Line2D([], [], marker="o", ls="", markersize=6, markerfacecolor="#bdbdbd", markeredgecolor="black")
    h_peak = plt.Line2D([], [], marker="D", ls="", markersize=9, markerfacecolor="#bdbdbd", markeredgecolor="black")
    ax.set_xticks(x, names, fontsize=TICK_SIZE, rotation=40, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Tajima's D (block median)", fontsize=LABEL_SIZE)
    ax.tick_params(axis="y", labelsize=TICK_SIZE)
    ax.spines[["top", "right"]].set_visible(False)
    # legend to the right of the axes
    ax.legend([bp1["boxes"][0], bp2["boxes"][0], h_other, h_peak],
              ["random blocks,\ngenome-wide", "random blocks,\ncentromere-matched", "other fog peaks", "chr 6 fog peak"],
              frameon=False, fontsize=TICK_SIZE - 3, loc="center left", bbox_to_anchor=(1.01, 0.5), labelspacing=1.0)
    fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    print(f"Saved: {OUT_FIG}")

    # ---- per-population KDE of genome-wide window Tajima's D, peaks below ----
    allw = np.concatenate([plot[n]["windows"] for n in names])
    xs = np.linspace(allw.min(), np.percentile(allw, 99.8), 400)
    fig, axs = plt.subplots(len(names), 1, figsize=(9, 2.1 * len(names)), sharex=True,
                            gridspec_kw=dict(hspace=0.25))
    for ax, n in zip(axs, names):
        d = gaussian_kde(plot[n]["windows"])(xs)
        ax.fill_between(xs, d, color=POP_COLOR[n], alpha=0.6, linewidth=0)
        ax.plot(xs, d, color="black", lw=1)
        top = d.max()
        o = plot[n]["other"]
        ax.scatter(o, -0.18 * top + rng.uniform(-0.05, 0.05, len(o)) * top, s=30, color=POP_COLOR[n],
                   edgecolor="black", linewidth=0.5, alpha=OTHER_ALPHA, zorder=3)
        ax.scatter(plot[n]["peak"], -0.18 * top, s=PEAK_SIZE * 0.65, marker="D", color=POP_COLOR[n], edgecolor="black",
                   linewidth=1.3, zorder=4)
        ax.set_ylim(-0.3 * top, top * 1.08)
        ax.set_yticks([])
        ax.text(0.99, 0.92, n, transform=ax.transAxes, fontsize=TICK_SIZE, va="top", ha="right")
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="x", labelsize=TICK_SIZE)
    axs[-1].set_xlabel("Tajima's D", fontsize=LABEL_SIZE)
    fig.text(0.07, 0.5, "Density (genome-wide windows)", rotation=90, va="center", ha="center", fontsize=LABEL_SIZE - 4)
    fig.savefig(OUT_FIG_KDE, dpi=200, bbox_inches="tight")
    print(f"Saved: {OUT_FIG_KDE}")


if __name__ == "__main__":
    main()
