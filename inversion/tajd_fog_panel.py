#!/usr/bin/env python3
"""
Stand-alone panel (alternative x-axis for the bottom-right panel of
DEINVI_0006_stacked.py): Delta Tajima's D (fog-peak windows = triangles,
other genome-wide fog outlier windows = circles; one point per locality,
coloured by locality) against each locality's long-term fog -- the
env_var.txt `fog` column (2000-2022 June-September mean of the MODIS-based
monthly fog rasters, env/fog/; Werner et al.), averaged over the plants in
that locality's Tajima's D sample (same bamlists as
plot_TajD_extreme_combined.R), untrustworthy-coordinate plants excluded.

Usage: tajd_fog_panel.py
Output: OUT
"""

import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma")
import temporal_null as tn                      # BAD_COORD_IDS

TAJD = "/private/groups/meyerlab/esteiger/selection/gemma/TajD_delta_vs_fogDetrend_fogPeak.tsv"
ENV = "/private/groups/meyerlab/esteiger/env_var.txt"
BAMLIST_DIR = "/private/groups/meyerlab/esteiger/mapping/merged"
TAJD_BAMLISTS = {"CDFW": "bamlist_4x_cdfw.txt", "Strauss": "bamlist_4x_strauss.txt",
                 "Point_Conception": "bamlist_4x_conception.txt", "VAFB_Lion_Head": "bamlist_4x_lionhead.txt",
                 "VAFB_Arguello": "bamlist_4x_arguello.txt"}
OUT = "/private/groups/meyerlab/esteiger/inversions/TajD_delta_vs_longterm_fog.jpg"
LABEL_SIZE, TICK_SIZE = 24, 19                  # as in DEINVI_0006_stacked.py
SHAPES = {"Outliers": "o", "Fog peak outliers": "^"}


def main():
    env = pd.read_csv(ENV, sep="\t", dtype={"SampleID": str}).drop_duplicates("SampleID").set_index("SampleID")
    pop_fog = {}
    for pop, f in TAJD_BAMLISTS.items():
        ids = [l.strip().replace("_rmdup.bam", "") for l in open(f"{BAMLIST_DIR}/{f}") if l.strip()]
        ids = [i for i in ids if i in env.index and i not in tn.BAD_COORD_IDS]
        pop_fog[pop] = env.loc[ids, "fog"].astype(float).mean()
    print("long-term fog (days/month) by locality: " + ", ".join(f"{k} {v:.2f}" for k, v in pop_fog.items()))
    tajd = pd.read_csv(TAJD, sep="\t")
    tajd["lt_fog"] = tajd.population.map(pop_fog)
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    ax.axhline(0, color="grey", linestyle="--", linewidth=1.5, zorder=0)
    for cat, mk in SHAPES.items():
        t = tajd[tajd.category == cat]
        ax.scatter(t.lt_fog, t.delta_TajD, marker=mk, s=230, c=t.color, edgecolor="black", linewidth=1.3,
                   alpha=0.9, zorder=3)
    ax.set_xlabel("Long-term fog\n(days/month)", fontsize=LABEL_SIZE)
    ax.set_ylabel("Δ Tajima's D", fontsize=LABEL_SIZE)
    ax.tick_params(labelsize=TICK_SIZE)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(tajd[["population", "category", "delta_TajD", "lt_fog"]].round(3).to_string(index=False))
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()
