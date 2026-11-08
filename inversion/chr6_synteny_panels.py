#!/usr/bin/env python3
"""
Supplementary figure for DEINVI_00006 (chr 6), x-aligned with
DEINVI_0006_stacked.py; grey box = centromere call (call_centromeres.py).

A: fog GWAS Manhattan (points exported by annotated_manhattan.R; dashed line =
   Bonferroni threshold).
B, C: winPCA PC1 (all-locality 5% MAF GL run, sign-flipped) and proportion of
   heterozygous sites (all-locality >= 8x GT run), one line per individual,
   coloured by fog; each locality capped at PER_SITE_CAP individuals for plotting.
D, E: GENESPACE syntenic hits (blkID assigned) with sunflower and lettuce; one
   band per partner chromosome with >= MIN_HITS hits.

Usage: chr6_synteny_panels.py
Output: OUT
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = "/private/groups/meyerlab/esteiger"

GS = f"{BASE_DIR}/inversions/synteny/genespace"
GEMMA_DIR = f"{BASE_DIR}/selection/gemma"
WINPCA_DIR = f"{BASE_DIR}/inversions/all_localities"
WINPCA_META = f"{BASE_DIR}/inversions/metadata.tsv"
CENTROMERES = f"{BASE_DIR}/inversions/centromere/centromeres.tsv"
FAI = f"{BASE_DIR}/inversions/filtered_scaffolds.fai"
OUT = f"{BASE_DIR}/inversions/DEINVI_0006_synteny.jpg"
CHROM = "DEINVI_00006"
PC1 = f"{WINPCA_DIR}/DEINVI_00006_GL_maf5_flipped.pc_1.tsv.gz"
HETP = f"{WINPCA_DIR}/DEINVI_00006_GT_all8x.hetp.tsv.gz"
MAN_LABELS = ["GNS1"]                                    # as in DEINVI_0006_stacked.py
LABEL_TEXT = {"GNS1": "GNS1, SNF1"}                      # display text (the peak spans both genes)
LABEL_SIZE, TICK_SIZE, PANEL_LETTER_SIZE = 24, 19, 28
SHADE_ALPHA, MIN_HITS, Y_PAD = 0.25, 20, 0.01
PER_SITE_CAP, SUBSAMPLE_SEED = 15, 1
LINE_WIDTH, LINE_ALPHA, CBAR_ALPHA = 1.8, 0.45, 0.6
# panel heights, in units of one partner-chromosome band
MAN_HEIGHT, PC1_HEIGHT, HET_HEIGHT, GAP_HEIGHT = 1.5, 2.0, 1.4, 0.35
LETTUCE = {f"CM0225{18 + i}.2": str(1 + i) for i in range(9)}   # CM022518.2 = Lsat chr 1, ...
SPECIES = [("Helianthus_annuus", "Sunflower"), ("Lactuca_sativa", "Lettuce")]
WINPCA_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "winpca_r", list(zip([0, 0.25, 0.5, 0.75, 1], ["#FFFF00", "#FF0000", "#800080", "#0000FF", "#000080"])))


def hits(sp):
    d = pd.read_csv(f"{GS}/syntenicHits/{sp}_vs_Deinandra.synHits.txt.gz", sep="\t", low_memory=False)
    d = d[(d.chr2 == CHROM) & d.blkID.notna()].copy()
    d["x"] = (d.start2 + d.end2) / 2e6
    d["y"] = (d.start1 + d.end1) / 2e6
    d["partner"] = d.chr1.astype(str).map(LETTUCE).fillna(d.chr1.astype(str))
    return d


def manhattan(ax):
    man = pd.read_csv(f"{GEMMA_DIR}/fog/manhattan_{CHROM}_fog.tsv.gz", sep="\t")
    lab = pd.read_csv(f"{GEMMA_DIR}/fog/manhattan_labels_{CHROM}_fog.tsv", sep="\t")
    lab = lab[lab["Label"].isin(MAN_LABELS)]
    man["logp"], lab["logp"] = -np.log10(man["p_final"]), -np.log10(lab["p_final"])
    ax.scatter(man["Position"] / 1e6, man["logp"], s=14, color="black", alpha=0.65,
               edgecolor="none", rasterized=True)
    ax.axhline(-np.log10(man["thresh_bonf"].iloc[0]), color="grey", linestyle="--", linewidth=1.5)
    for _, r in lab.iterrows():
        ax.annotate(LABEL_TEXT.get(r["Label"], r["Label"]), xy=(r["Position"] / 1e6, r["logp"]), xytext=(3, 6),
                    textcoords="offset points", ha="left", va="center", fontsize=LABEL_SIZE,
                    fontstyle="italic", arrowprops=dict(arrowstyle="-", color="grey", linewidth=0.8))
    ax.set_ylim(0, man["logp"].max() * (1 + Y_PAD))
    ax.set_ylabel("-log10\np-val", fontsize=LABEL_SIZE)


def winpca_lines(pc1, hetp):
    """Per-site subsample, fog values, colour normalization and drawing order."""
    meta = pd.read_csv(WINPCA_META, sep="\t", index_col=0)
    samples = pc1.columns.union(hetp.columns)
    fog_val = -pd.to_numeric(meta["fog"], errors="coerce").reindex(samples)
    site = meta["Clusters"].reindex(samples)
    rng = np.random.default_rng(SUBSAMPLE_SEED)
    in_both = set(pc1.columns) & set(hetp.columns)
    drawn = set()
    for _, members in site.groupby(site):
        both = sorted(s for s in members.index if s in in_both)
        rest = sorted(s for s in members.index if s not in in_both)
        picked = list(rng.permutation(both))[:PER_SITE_CAP]
        if len(picked) < PER_SITE_CAP:
            picked += list(rng.permutation(rest))[:PER_SITE_CAP - len(picked)]
        drawn.update(picked)
    norm = matplotlib.colors.Normalize(vmin=fog_val.min(), vmax=fog_val.max())
    keep = fog_val[fog_val.index.isin(drawn)]
    order = list(keep.index[keep.isna()]) + list(keep.dropna().sort_values(ascending=False).index)
    return fog_val, norm, order


def line_panel(ax, df, ylabel, fog_val, norm, order):
    x = df.index.values / 1e6
    for s in [s for s in order if s in df.columns]:
        v = fog_val[s]
        ax.plot(x, df[s].values, color="lightgrey" if pd.isna(v) else WINPCA_CMAP(norm(v)),
                linewidth=LINE_WIDTH, alpha=LINE_ALPHA, rasterized=True)
    ax.margins(y=Y_PAD)
    ax.set_ylabel(ylabel, fontsize=LABEL_SIZE)


def main():
    data = {sp: hits(sp) for sp, _ in SPECIES}
    xmax = pd.read_csv(FAI, sep="\t", header=None).set_index(0).loc[CHROM, 1] / 1e6
    cen = pd.read_csv(CENTROMERES, sep="\t").query("confident").set_index("chr").loc[CHROM]
    shade = [(cen.cen_start / 1e6, cen.cen_end / 1e6)]
    pc1 = pd.read_csv(PC1, sep="\t", index_col=0)
    hetp = pd.read_csv(HETP, sep="\t", index_col=0)
    fog_val, norm, order = winpca_lines(pc1, hetp)

    # one band per partner chromosome, ordered by chr 6 position of its first hits
    bands = {}
    for sp, _ in SPECIES:
        d = data[sp]
        n = d.groupby("partner").size()
        keep = n[n >= MIN_HITS].index
        bands[sp] = sorted(keep, key=lambda c: d.loc[d.partner == c, "x"].quantile(0.05))
        print(f"{sp}: bands {[(c, int(n[c])) for c in bands[sp]]}; "
              f"dropped (< {MIN_HITS} hits): {dict(n[n < MIN_HITS])}")

    ratios = [MAN_HEIGHT, PC1_HEIGHT, HET_HEIGHT, GAP_HEIGHT]
    for i, (sp, _) in enumerate(SPECIES):
        ratios += [1] * len(bands[sp]) + ([GAP_HEIGHT] if i < len(SPECIES) - 1 else [])
    fig = plt.figure(figsize=(14, 1.55 * sum(ratios) + 1.4))
    gs = fig.add_gridspec(len(ratios), 1, height_ratios=ratios, hspace=0.12)

    man_ax = fig.add_subplot(gs[0])
    manhattan(man_ax)
    pc1_ax = fig.add_subplot(gs[1], sharex=man_ax)
    line_panel(pc1_ax, pc1, "PC1", fog_val, norm, order)
    het_ax = fig.add_subplot(gs[2], sharex=man_ax)
    line_panel(het_ax, hetp, "Proportion\nhet. sites", fog_val, norm, order)
    top_axes = [man_ax, pc1_ax, het_ax]
    all_axes, groups, row = list(top_axes), [], 4
    for si, (sp, name) in enumerate(SPECIES):
        d, grp = data[sp], []
        for c in bands[sp]:
            ax = fig.add_subplot(gs[row], sharex=man_ax); row += 1
            t = d[d.partner == c]
            ax.scatter(t.x, t.y, s=22, c="black", edgecolor="none", alpha=0.9, zorder=3)
            lo, hi = t.y.min(), t.y.max(); pad = 0.08 * (hi - lo + 1)
            ax.set_ylim(lo - pad, hi + pad)
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
            ax.text(1.01, 0.5, f"chr {c}", transform=ax.transAxes, fontsize=TICK_SIZE,
                    va="center", ha="left")
            grp.append(ax); all_axes.append(ax)
        groups.append(grp)
        row += 1
    for ax in all_axes:
        for a, b in shade:
            ax.axvspan(a, b, color="#808080", alpha=SHADE_ALPHA, lw=0, zorder=0)
        ax.tick_params(axis="both", labelsize=TICK_SIZE)
        ax.tick_params(axis="x", labelbottom=ax is all_axes[-1])
        ax.spines[["top", "right"]].set_visible(False)
    man_ax.set_xlim(0, xmax)
    all_axes[-1].set_xlabel(f"Position on chromosome {CHROM} (Mb)", fontsize=LABEL_SIZE)
    fig.subplots_adjust(left=0.12, right=0.9, top=0.98, bottom=0.06)

    # y-titles and panel letters
    for ax in top_axes:
        ax.yaxis.set_label_coords(-0.06, 0.5)
    het_ax.yaxis.set_label_coords(-0.048, 0.5)    # two-line title: keep clear of the panel letter
    for (sp, name), grp in zip(SPECIES, groups):
        top, bot = grp[0].get_position().y1, grp[-1].get_position().y0
        fig.text(0.05, (top + bot) / 2, f"{name} (Mb)", rotation=90, va="center", ha="center",
                 fontsize=LABEL_SIZE)
    for letter, ax in zip("ABCDE", top_axes + [g[0] for g in groups]):
        fig.text(0.005, ax.get_position().y1, letter, fontsize=PANEL_LETTER_SIZE, va="top", ha="left")
    # fog colorbar spanning B-C
    tp, bp = pc1_ax.get_position(), het_ax.get_position()
    cax = fig.add_axes([tp.x1 + 0.01, bp.y0 + 0.02 * bp.height, 0.012, tp.y1 - bp.y0 - 0.04 * bp.height])
    cbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=WINPCA_CMAP), cax=cax, alpha=CBAR_ALPHA)
    cbar.set_label("Fog (days per month)", fontsize=LABEL_SIZE)
    cbar.ax.tick_params(labelsize=TICK_SIZE)
    fig.savefig(OUT, dpi=200)
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()
