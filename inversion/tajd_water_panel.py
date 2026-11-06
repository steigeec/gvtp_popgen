#!/usr/bin/env python3
"""
Stand-alone two-panel figure (local-adaptation framing, alternative to the
bottom-right panel of DEINVI_0006_stacked.py): Delta Tajima's D (fog-peak
windows = triangles, other genome-wide fog outlier windows = circles; one
point per locality, coloured by locality) against water availability.

Hypothesis: selection for fog-favourable alleles should be strongest where
fog is plentiful relative to other water. Two locality-level axes, averaged
over the plants in each locality's Tajima's D sample (bamlists as in
plot_TajD_extreme_combined.R; untrustworthy-coordinate plants excluded):
  A. fog relative to rainfall: long-term fog (env_var `fog`, 2000-2022
     June-September mean, days/month) per 100 mm of mean annual
     precipitation (env_var `BIO12`).
  B. climatic water deficit (env_var `BCM_30y`, 30-year mean, mm): overall
     water stress, independent of fog.
Spearman correlations are printed and shown, but with five localities they are
descriptive (|rho| = 1 is needed for p < 0.05, two-sided).

Usage: tajd_water_panel.py
Output: OUT
"""

import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sys.path.insert(0, "/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma")
import temporal_null as tn                      # BAD_COORD_IDS

TAJD = "/private/groups/meyerlab/esteiger/selection/gemma/TajD_delta_vs_fogDetrend_fogPeak.tsv"
ENV = "/private/groups/meyerlab/esteiger/env_var.txt"
BAMLIST_DIR = "/private/groups/meyerlab/esteiger/mapping/merged"
TAJD_BAMLISTS = {"CDFW": "bamlist_4x_cdfw.txt", "Strauss": "bamlist_4x_strauss.txt",
                 "Point_Conception": "bamlist_4x_conception.txt", "VAFB_Lion_Head": "bamlist_4x_lionhead.txt",
                 "VAFB_Arguello": "bamlist_4x_arguello.txt"}
OUT = "/private/groups/meyerlab/esteiger/inversions/TajD_delta_vs_water.jpg"
LABEL_SIZE, TICK_SIZE = 24, 19                  # as in DEINVI_0006_stacked.py
SHAPES = {"Outliers": "o", "Fog peak outliers": "^"}


def main():
    env = pd.read_csv(ENV, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    loc = {}
    for pop, f in TAJD_BAMLISTS.items():
        ids = [l.strip().replace("_rmdup.bam", "") for l in open(f"{BAMLIST_DIR}/{f}") if l.strip()]
        ids = [i for i in ids if i in env.index and i not in tn.BAD_COORD_IDS]
        e = env.loc[ids, ["fog", "BIO12", "BCM_30y"]].apply(pd.to_numeric, errors="coerce").mean()
        loc[pop] = dict(fog_per_100mm=100 * e.fog / e.BIO12, cwd=e.BCM_30y, fog=e.fog, precip=e.BIO12)
    loc = pd.DataFrame(loc).T
    print(loc.round(3).to_string())
    tajd = pd.read_csv(TAJD, sep="\t").join(loc, on="population")

    axes_def = [("fog_per_100mm", "Fog relative to rainfall\n(fog days per 100 mm)"),
                ("cwd", "Climatic water deficit\n(mm)")]
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.2), gridspec_kw=dict(wspace=0.35))
    for ax, (col, xl) in zip(axs, axes_def):
        ax.axhline(0, color="grey", linestyle="--", linewidth=1.5, zorder=0)
        txt = []
        for cat, mk in SHAPES.items():
            t = tajd[tajd.category == cat]
            ax.scatter(t[col], t.delta_TajD, marker=mk, s=230, c=t.color, edgecolor="black", linewidth=1.3,
                       alpha=0.9, zorder=3)
            rho, p = spearmanr(t[col], t.delta_TajD)
            txt.append(f"{'peak' if cat.startswith('Fog') else 'other outliers'}: rho = {rho:+.2f}")
            print(f"{col:14s} {cat:18s} Spearman rho {rho:+.2f} (p = {p:.2f}, n = {len(t)})")
        ax.set_xlabel(xl, fontsize=LABEL_SIZE)
        ax.set_ylabel("Δ Tajima's D", fontsize=LABEL_SIZE)
        ax.tick_params(labelsize=TICK_SIZE)
        ax.spines[["top", "right"]].set_visible(False)
    for ax, l in zip(axs, "AB"):
        ax.text(-0.25, 1.02, l, transform=ax.transAxes, fontsize=28, va="bottom")
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()
