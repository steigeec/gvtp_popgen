#!/usr/bin/env python
"""
Fit FEEMS for one sampling year. Lambda is chosen by leave-one-node-out
cross-validation unless given. The edge colour map passes through grey
instead of white so edges near the average rate stay visible.

Usage: python FEEMS_fit.py <plink_prefix> [lambda]
Output: <plink_prefix>_lambda<lambda>.jpg/.pdf; with CV also <plink_prefix>_cv.tsv and <plink_prefix>_cv.jpg
"""

import sys

import numpy as np
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import cartopy.crs as ccrs
from feems import SpatialGraph, Viz
from feems.cross_validation import run_cv

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"

EDGE_WIDTH = 2.0
LAMB_GRID = np.geomspace(1e-6, 1e4, 120)[::-1]     # large to small for warm starts
OBS_NODE_SIZE = 9
SAMPLE_PT_SIZE = 20
# orange (low) to grey to blue (high)
EDGE_COLORS = ["#6b2600", "#a33d00", "#d9620b", "#f29a4a", "#a8a8a8",
               "#5bb8d6", "#1f8fb8", "#006e94", "#004a66"]

if len(sys.argv) < 2:
    print("Usage: python FEEMS_fit.py <plink_prefix> [lambda]")
    sys.exit(1)
prefix = sys.argv[1]
lamb = float(sys.argv[2]) if len(sys.argv) > 2 else None

# mean-impute genotypes, keep variable sites
bim, fam, G = read_plink(prefix)
genotypes = SimpleImputer(missing_values=np.nan, strategy="mean").fit_transform(np.array(G).T)
genotypes = genotypes[:, np.var(genotypes, axis=0) > 0]
print(f"{prefix}: n_samples={genotypes.shape[0]}, n_snps={genotypes.shape[1]}", flush=True)

coord = np.loadtxt(f"{prefix}.coord")
outer = np.loadtxt("outer.csv")
grid = np.loadtxt("grid.csv")
edges = np.loadtxt("edges.csv", dtype=int)
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

if lamb is None:
    cv_err = run_cv(sp_graph, LAMB_GRID, n_folds=sp_graph.n_observed_nodes, factr=1e10)
    mean_err = np.mean(cv_err, axis=0)[:, 0]
    lamb = float(LAMB_GRID[np.argmin(mean_err)])
    np.savetxt(f"{prefix}_cv.tsv", np.column_stack([LAMB_GRID, mean_err]), delimiter="\t",
               header="lambda\tmean_cv_l2_error", comments="")
    fig, ax = plt.subplots(dpi=300)
    ax.plot(np.log10(LAMB_GRID), mean_err, ".", color="black")
    ax.axvline(np.log10(lamb), color="#d9620b")
    ax.set_xlabel("log10(lambda)"); ax.set_ylabel("Mean L2 CV error")
    ax.set_title(f"{prefix}: lambda = {lamb:.4g}", fontsize=9)
    fig.savefig(f"{prefix}_cv.jpg", bbox_inches="tight"); plt.close(fig)
    print(f"cross-validated lambda = {lamb:.6g} ({sp_graph.n_observed_nodes} folds)", flush=True)
sp_graph.fit(lamb)

projection = ccrs.EquidistantConic(central_longitude=coord[:, 0].mean(),
                                   central_latitude=coord[:, 1].mean())
fig = plt.figure(figsize=(8, 6), dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=EDGE_WIDTH,
        edge_alpha=1, edge_zorder=100, sample_pt_size=SAMPLE_PT_SIZE,
        obs_node_size=OBS_NODE_SIZE, sample_pt_color="black", cbar_font_size=10)
# palette swapped on the Viz object so the installed feems package is untouched
v.edge_cmap = clr.LinearSegmentedColormap.from_list("feems_no_white", EDGE_COLORS, N=256)
v.draw_map()
v.draw_samples()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
v.draw_edge_colorbar()
ax.set_aspect("equal")
tag = f"{lamb:.4g}"
for ext in ("jpg", "pdf"):
    plt.savefig(f"{prefix}_lambda{tag}.{ext}", bbox_inches="tight")
print(f"Saved: {prefix}_lambda{tag}.jpg / .pdf")
