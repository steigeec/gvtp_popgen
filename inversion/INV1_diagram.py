#!/usr/bin/env python3
"""
Evidence map of DEINVI_00006 ~10-52 Mb (the region previously called
"INV1"/"INV2"). Only features supported by this project's analyses (Oct 2026)
are drawn; neither "INV1" nor "INV2" is drawn, since no inversion is
read-confirmed. Top: whole region; bottom: zoom on the fog outlier peak.

Tracks and sources
  Synteny (GENESPACE vs H. annuus, H. anomalus, Lactuca sativa):
      0.3-31.9 Mb collinear with one ancestral chromosome segment;
      19.34-20.83 Mb reversed relative to all three; 32-48.7 Mb no syntenic
      genes; >= 48.76 Mb a different ancestral segment.
  LINE content: 45.5-48 Mb enriched (~5-12% vs ~1% genome-wide; RepeatMasker).
  Population haplotype structure (winPCA 5% MAF, cluster consistency):
      elevated PC1 % variance explained ~16-49.3 Mb; blocks where samples
      keep a stable three-group split (core-window extents, all-locality run);
      sharp boundary ~49.1-49.4 Mb (PCA, heterozygosity, gene density).
  Read-level junctions (split reads / read pairs):
      polymorphic -- present in one PC group, absent in the other (population
        short reads, tight coordinates; structure unresolved);
      reference discrepancy -- present in most individuals and in the
        reference plant's own long reads, i.e. the assembly order/orientation
        likely differs from the population there.
  COPG2 rearrangement (48.8845-48.909 Mb): read-confirmed, genotyped,
      Mendelian; Strauss frequency 0.06, ~fixed at Gaviota.
  Fog outlier peak (49.031-49.146 Mb; GNS1, SNF1) and polymorphic
      structural variants found by peak_sv_scan.py (all samples, unbiased).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Arc, Rectangle

OUT = "/private/groups/meyerlab/esteiger/inversions/INV1_diagram.jpg"
SV_EVENTS = "/private/groups/meyerlab/esteiger/inversions/breakpoints/peak_sv_events.tsv"
FS = 12
C = {"syn": "#9e9e9e", "syn2": "#5f5f5f", "nosyn": "#ffffff", "line": "#7b3294", "pc": "#fdd49e",
     "blk": "#d95f02", "poly": "#1b9e77", "asm": "#999999", "copg2": "#000000", "peak": "#f0c400", "sv": "#e7298a"}

SYNTENY = [(10.0, 31.9, "syn"), (48.76, 52.0, "syn2")]
REVERSED = (19.34, 20.83)
LINE_ENRICHED = (45.5, 48.0)
PC1_ELEVATED = (16.0, 49.3)
BLOCKS = [(16.0, 18.1), (19.1, 21.4), (21.4, 23.1), (24.3, 27.3), (27.7, 29.0), (31.9, 36.6),
          (44.2, 45.0), (45.5, 46.5), (47.5, 49.1)]
BOUNDARY = (49.1, 49.4)
POLY_JUNCTIONS = [(16.66, 18.83), (19.80, 37.48), (18.79, 37.99)]
ASM_JUNCTIONS = [(16.03, 16.68), (21.00, 45.30), (22.14, 24.90), (31.35, 37.47)]
COPG2 = (48.8845, 48.909)
PEAK = (49.031, 49.146)
GENES = {"GNS1": (49.0998, 49.1045), "SNF1": (49.1214, 49.1240), "COPG2": (48.8777, 48.8881)}


def peak_svs():
    """Polymorphic SV events from the unbiased scan (peak_sv_scan.py, all 339
    samples; clip positions within 100 bp merged into one event): carried by
    >= 5 individuals and < 95% of genotyped individuals, genotypes consistent
    with Hardy-Weinberg in Strauss."""
    try:
        e = pd.read_csv(SV_EVENTS, sep="\t")
    except FileNotFoundError:
        return pd.DataFrame(columns=["start", "kind", "rho_fav"])
    return e[e.cls == "polymorphic (HWE ok)"]


def track_label(ax, y, text):
    ax.text(-0.005, y, text, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=FS)


fig, (ax, az) = plt.subplots(2, 1, figsize=(15, 11.5), gridspec_kw={"height_ratios": [1.3, 1], "hspace": 0.32})

# ── whole region ─────────────────────────────────────────────────────────────
X0, X1 = 10, 52
ax.set_xlim(X0, X1); ax.set_ylim(-0.4, 8.2)
ax.set_yticks([]); ax.spines[["left", "right", "top"]].set_visible(False)
ax.set_xticks(range(10, 53, 2)); ax.tick_params(axis="x", labelsize=FS)
ax.set_xlabel("Position on chromosome 6 (Mb)", fontsize=FS + 1)

y = 7.4; track_label(ax, y, "Synteny with\nsunflower / lettuce")
ax.add_patch(Rectangle((X0, y - 0.25), X1 - X0, 0.5, facecolor=C["nosyn"], edgecolor="black", lw=0.6))
for a, b, k in SYNTENY:
    ax.add_patch(Rectangle((a, y - 0.25), b - a, 0.5, color=C[k], lw=0))
ax.add_patch(Rectangle((REVERSED[0], y - 0.25), REVERSED[1] - REVERSED[0], 0.5, facecolor="white", edgecolor="black", hatch="////", lw=0.8))
ax.text(20.1, y + 0.35, "reversed vs relatives", ha="center", va="bottom", fontsize=FS - 2)
ax.text(40.3, y, "no syntenic genes", ha="center", va="center", fontsize=FS - 2)
ax.text(21, y - 0.42, "ancestral segment 1", ha="center", va="top", fontsize=FS - 2, color="#444444")
ax.text(50.4, y - 0.42, "segment 2", ha="center", va="top", fontsize=FS - 2, color="#444444")

y = 6.0; track_label(ax, y, "LINE-enriched")
ax.add_patch(Rectangle((LINE_ENRICHED[0], y - 0.2), LINE_ENRICHED[1] - LINE_ENRICHED[0], 0.4, color=C["line"], lw=0))

y = 4.8; track_label(ax, y, "Population haplotype\nstructure (winPCA)")
ax.add_patch(Rectangle((PC1_ELEVATED[0], y - 0.35), PC1_ELEVATED[1] - PC1_ELEVATED[0], 0.7, color=C["pc"], lw=0))
for a, b in BLOCKS:
    ax.add_patch(Rectangle((a, y - 0.22), b - a, 0.44, facecolor=C["blk"], lw=0))
ax.add_patch(Rectangle((BOUNDARY[0], y - 0.5), BOUNDARY[1] - BOUNDARY[0], 1.0, facecolor="none", edgecolor="black", lw=2))
ax.text(49.25, y + 0.6, "sharp\nboundary", ha="center", va="bottom", fontsize=FS - 2)

y = 2.7; track_label(ax, y, "Read-level\njunctions")
ax.plot([X0, X1], [y, y], color="#dddddd", lw=1)
for (a, b), col, ls in [(j, C["poly"], "-") for j in POLY_JUNCTIONS] + [(j, C["asm"], "--") for j in ASM_JUNCTIONS]:
    h = 0.35 + 0.045 * (b - a)
    ax.add_patch(Arc(((a + b) / 2, y), b - a, 2 * h, theta1=0, theta2=180, color=col, lw=2, ls=ls))
ax.plot([sum(COPG2) / 2], [y], marker="v", color=C["copg2"], ms=11)
ax.text(sum(COPG2) / 2, y - 0.3, "COPG2\nrearr.", ha="center", va="top", fontsize=FS - 2)

y = 0.6; track_label(ax, y, "fog GWAS")
ax.add_patch(Rectangle((PEAK[0], y - 0.3), PEAK[1] - PEAK[0] + 0.15, 0.6, color=C["peak"], lw=0))
ax.text(PEAK[1] + 0.35, y, "outlier peak", ha="left", va="center", fontsize=FS - 1)
ax.add_patch(Rectangle((48.8, -0.35), 0.5, 8.3, facecolor="none", edgecolor="#666666", lw=1, ls=":"))

handles = [Rectangle((0, 0), 1, 1, color=C["syn"]), Rectangle((0, 0), 1, 1, color=C["pc"]),
           Rectangle((0, 0), 1, 1, color=C["blk"]), plt.Line2D([], [], color=C["poly"], lw=2),
           plt.Line2D([], [], color=C["asm"], lw=2, ls="--")]
labels = ["collinear with relatives", "elevated PC1 variance explained", "stable three-group haplotype block",
          "polymorphic junction (short reads)", "reference discrepancy (most individuals + ref. long reads)"]
ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=3, frameon=False, fontsize=FS - 1)

# ── zoom: fog peak ───────────────────────────────────────────────────────────
Z0, Z1 = 48.8, 49.3
az.set_xlim(Z0, Z1); az.set_ylim(-0.05, 4.6)
az.set_yticks([]); az.spines[["left", "right", "top"]].set_visible(False)
az.tick_params(axis="x", labelsize=FS)
az.set_xlabel("Position on chromosome 6 (Mb) - zoom on the fog outlier peak", fontsize=FS + 1)

y = 3.9; track_label(az, y, "Genes")
for g, (a, b) in GENES.items():
    az.add_patch(Rectangle((a, y - 0.15), b - a, 0.3, color="black", lw=0))
    az.text((a + b) / 2, y + 0.25, g, ha="center", va="bottom", fontsize=FS - 1, style="italic")

y = 2.9; track_label(az, y, "fog GWAS")
az.add_patch(Rectangle((PEAK[0], y - 0.25), PEAK[1] - PEAK[0], 0.5, color=C["peak"], lw=0))
az.text(sum(PEAK) / 2, y, "outlier peak (19 SNPs)", ha="center", va="center", fontsize=FS - 1)
az.add_patch(Rectangle((BOUNDARY[0], y - 0.45), BOUNDARY[1] - BOUNDARY[0], 0.9, facecolor="none", edgecolor="black", lw=2))
az.text(49.2, y + 0.55, "sharp haplotype-structure boundary", ha="center", va="bottom", fontsize=FS - 2)

y = 1.6; track_label(az, y, "Polymorphic\nstructural variants\n(short reads, all samples)")
az.add_patch(Rectangle((COPG2[0], y - 0.22), COPG2[1] - COPG2[0], 0.44, color=C["copg2"], lw=0))
az.text(sum(COPG2) / 2, y - 0.35, "COPG2 rearrangement\n(Strauss 0.06; ~fixed Gaviota)", ha="center", va="top", fontsize=FS - 2)
svs = peak_svs()
kcol = {"transposon-like insertion": "#e7298a", "local rearrangement": "#1b9e77", "unaligned clip": "#7570b3"}
mk = {"transposon-like insertion": "^", "local rearrangement": "s", "unaligned clip": "o"}
y_rug = 0.55
track_label(az, y_rug, "all polymorphic\nSV events (n=%d)" % len(svs))
az.add_patch(Rectangle((48.95, y_rug - 0.32), 0.30, 0.64, facecolor="#f2f2f2", edgecolor="#999999", lw=0.8, zorder=0))
az.text(48.95, y_rug - 0.36, "scanned interval (48.95-49.25 Mb)", ha="left", va="top", fontsize=FS - 3, color="#666666")
for r in svs.itertuples():
    az.plot([r.start / 1e6] * 2, [y_rug - 0.25, y_rug + 0.25], color=kcol.get(r.kind, "grey"), lw=0.8)
linked = svs[svs.rho_fav.abs() >= 0.5]
for r in linked.itertuples():
    az.plot([r.start / 1e6], [y], marker=mk.get(r.kind, "o"), color=kcol.get(r.kind, "grey"), ms=11, mec="black", mew=0.7, ls="none")
az.text(Z1, y + 0.45, f"large symbols: the {len(linked)} events linked to the fog-peak haplotype (|rho| >= 0.5)",
        ha="right", va="bottom", fontsize=FS - 2, color="#444444")
az.legend([plt.Line2D([], [], marker=mk[k], color=kcol[k], mec="black", ls="none", ms=9) for k in mk],
          list(mk.keys()), loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False, fontsize=FS - 1)

plt.savefig(OUT, dpi=200, bbox_inches="tight", pad_inches=0.15)
print(f"Saved: {OUT} ({len(svs)} polymorphic peak-region SVs drawn)")
