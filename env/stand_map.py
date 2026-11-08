#!/usr/bin/env python3
"""
Map of sampling stands (field_data.txt OHANA_order), one panel per genetic
cluster, with plants coloured and labelled by stand. Tranquillon Mtn
(Strauss) gets a larger panel. The 13 individuals with unreliable
coordinates (BAD_COORD_IDS) are left out.

Usage: stand_map.py
Output: OUT
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

BASE_DIR = "/private/groups/meyerlab/esteiger"

FIELD = f"{BASE_DIR}/field_data.txt"
OUT = f"{BASE_DIR}/env/stand_map.jpg"
BAD_COORD_IDS = {"4548", "4549", "4550", "4551", "4552", "4553", "4554",
                 "4555", "4556", "4557", "4559", "4561", "4562"}
NAMES = {"Strauss": "Tranquillon Mtn", "CDFW": "Gaviota", "Dangermond": "Dangermond",
         "Point_Conception": "Pt. Conception", "VAFB_Arguello": "Pt. Arguello", "VAFB_Lion_Head": "Lion's Head"}
BIG = "Strauss"
PAD_FRAC, SCALE_M = 0.08, 100


def draw(ax, d, title, point_size, label_size):
    stands = sorted(d.OHANA_order.unique())
    cmap = plt.get_cmap("tab20" if len(stands) > 10 else "tab10")
    # equal metres in x and y; square window centred on the full extent, with padding
    lat0 = (d.Lat.min() + d.Lat.max()) / 2; lon0 = (d.Long.min() + d.Long.max()) / 2
    k = np.cos(np.radians(lat0))
    xr, yr = np.ptp(d.Long) * k, np.ptp(d.Lat)
    half = max(xr, yr, 0.0005) / 2 * (1 + 2 * PAD_FRAC)
    ax.set_xlim(lon0 - half / k, lon0 + half / k)
    ax.set_ylim(lat0 - half, lat0 + half)
    for i, s in enumerate(stands):
        g = d[d.OHANA_order == s]
        col = cmap(i % cmap.N)
        ax.scatter(g.Long, g.Lat, s=point_size, color=col, edgecolor="black", linewidth=0.3, zorder=3)
        # label above the northernmost plant so it does not cover points
        ax.text(g.Long.median(), g.Lat.max() + 0.02 * half, str(s), color=col, fontsize=label_size,
                fontweight="bold", ha="center", va="bottom", zorder=4,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")])
    ax.set_aspect(1 / k)
    # scale bar
    span_m = 2 * half * 110570
    bar = SCALE_M if span_m < 3000 else 500
    x0 = ax.get_xlim()[0] + 0.05 * np.diff(ax.get_xlim())[0]; y0 = ax.get_ylim()[0] + 0.05 * np.diff(ax.get_ylim())[0]
    ax.plot([x0, x0 + bar / (111320 * k)], [y0, y0], color="black", linewidth=2.5)
    ax.text(x0 + bar / (111320 * k) / 2, y0 + 0.02 * np.diff(ax.get_ylim())[0], f"{bar} m", ha="center",
            va="bottom", fontsize=9)
    ax.set_title(f"{title} ({len(d)} plants, {len(stands)} stands)", fontsize=12)
    ax.ticklabel_format(useOffset=False)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.tick_params(labelsize=7)
    ax.set_xlabel("Longitude", fontsize=8); ax.set_ylabel("Latitude", fontsize=8)


def main():
    f = pd.read_csv(FIELD, sep="\t", dtype={"SampleID": str})
    f = f[~f.SampleID.isin(BAD_COORD_IDS)].dropna(subset=["Lat", "Long", "OHANA_order"])
    f["OHANA_order"] = f.OHANA_order.astype(int)
    others = [c for c in sorted(f.Clusters.unique()) if c != BIG]
    fig = plt.figure(figsize=(20, 10))
    gs = fig.add_gridspec(2, 5, wspace=0.35, hspace=0.3)
    draw(fig.add_subplot(gs[:, :2]), f[f.Clusters == BIG], NAMES.get(BIG, BIG), 14, 13)
    slots = [gs[0, 2], gs[0, 3], gs[0, 4], gs[1, 2], gs[1, 3], gs[1, 4]]
    for c, slot in zip(others, slots):
        draw(fig.add_subplot(slot), f[f.Clusters == c], NAMES.get(c, c), 22, 11)
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()
