#!/usr/bin/env python
"""
Fit FRAME (directed FEEMS; https://github.com/NovembreLab/frame) on the FEEMS
grid. Inputs: <prefix>.bed/.bim/.fam, <prefix>.coord, grid.csv, edges.csv.
lambda_m is chosen by 10-fold cross-validation unless given.

FRAME's rate m[i->j] is the backward-in-time lineage rate, so forward gene flow
is j to i. The edge table records both directions.

Usage: python FRAME_fit.py <plink_prefix> [lambda_m]   (no lambda = cross-validate)
Output: <prefix>_frame_cv.tsv/.jpg (with CV), <prefix>_frame_lambda<l>_edges.tsv,
        <prefix>_frame_lambda<l>_nodes.tsv, <prefix>_frame_lambda<l>.jpg/.pdf
"""

import os
import sys

os.environ["DISABLE_PANDERA_IMPORT_WARNING"] = "True"

import numpy as np
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from frame import SpatialDiGraph, Vis
from frame.cross_validation import run_cv

LAMB_WARMUP = 1e3
LAMB_GRID = np.geomspace(1e-3, 1e3, 20)[::-1]   # FRAME example workflow's grid
N_FOLDS = 10
RANDOM_STATE = 200

if len(sys.argv) < 2:
    print("Usage: python FRAME_fit.py <plink_prefix> [lambda_m]")
    sys.exit(1)
prefix = sys.argv[1]
lamb = float(sys.argv[2]) if len(sys.argv) > 2 else None

# mean-impute genotypes, keep variable sites
bim, fam, G = read_plink(prefix)
genotypes = SimpleImputer(missing_values=np.nan, strategy="mean").fit_transform(np.array(G).T)
genotypes = genotypes[:, np.var(genotypes, axis=0) > 0]
print(f"{prefix}: n_samples={genotypes.shape[0]}, n_snps={genotypes.shape[1]}", flush=True)

coord = np.loadtxt(f"{prefix}.coord")
grid = np.loadtxt("grid.csv")
# edges.csv is 1-based (FEEMS convention) but FRAME needs 0-based rows of grid.csv
edges = np.loadtxt("edges.csv", dtype=int) - 1
assert edges.min() == 0 and edges.max() == len(grid) - 1, "edges must index grid.csv rows (0-based)"
sp_digraph = SpatialDiGraph(genotypes, coord, grid, edges)

if lamb is None:
    cv_errs, _ = run_cv(sp_digraph, lamb_m_grid=LAMB_GRID, lamb_m_warmup=LAMB_WARMUP,
                        n_folds=N_FOLDS, factr=1e10, random_state=RANDOM_STATE)
    cv_errs = np.asarray(cv_errs).ravel()
    lamb = float(f"{LAMB_GRID[np.argmin(cv_errs)]:.3g}")
    np.savetxt(f"{prefix}_frame_cv.tsv", np.column_stack([LAMB_GRID, cv_errs]), delimiter="\t",
               header="lambda_m\tcv_error", comments="")
    fig, ax = plt.subplots(dpi=300)
    ax.plot(np.log10(LAMB_GRID), cv_errs, "o", color="black")
    ax.axvline(np.log10(lamb), color="#d9620b")
    ax.set_xlabel("log10(lambda_m)"); ax.set_ylabel("CV error")
    ax.set_title(f"{prefix}: FRAME lambda_m = {lamb:.3g}", fontsize=9)
    fig.savefig(f"{prefix}_frame_cv.jpg", bbox_inches="tight"); plt.close(fig)
    print(f"cross-validated lambda_m = {lamb:.3g}", flush=True)

# warm start from a heavily smoothed fit
sp_digraph.fit(lamb_m=LAMB_WARMUP, factr=1e10)
init = dict(logm_init=np.log(sp_digraph.m), logc_init=np.log(sp_digraph.c),
            trans_alpha_init=-np.log((1 / sp_digraph.alpha) - 1))
sp_digraph.fit(lamb_m=lamb, factr=1e7, **init)

# ---------- directed edge table ----------
# rows of M are lineage sources, columns lineage destinations
M = sp_digraph.M.tocoo()
tag = f"{lamb:.3g}"   # lambda in filenames so runs don't overwrite
with open(f"{prefix}_frame_lambda{tag}_edges.tsv", "w") as fh:
    fh.write("lineage_from\tlineage_to\tgeneflow_from\tgeneflow_to\t"
             "geneflow_from_lon\tgeneflow_from_lat\tgeneflow_to_lon\tgeneflow_to_lat\tm\n")
    for i, j, m in zip(M.row, M.col, M.data):
        fh.write(f"{i}\t{j}\t{j}\t{i}\t{grid[j, 0]:.6f}\t{grid[j, 1]:.6f}\t"
                 f"{grid[i, 0]:.6f}\t{grid[i, 1]:.6f}\t{m:.6g}\n")

n_hap = np.array([sp_digraph.nodes[k]["n_haps"] for k in range(sp_digraph.number_of_nodes())])
np.savetxt(f"{prefix}_frame_lambda{tag}_nodes.tsv",
           np.column_stack([np.arange(len(grid)), grid, n_hap, np.ravel(sp_digraph.gamma), np.ravel(sp_digraph.pi)]),
           delimiter="\t", header="node\tlon\tlat\tn_haps\tcoalescence_rate\tstationary_pi", comments="", fmt="%g")
print(f"lambda_m = {lamb:.3g}: {len(M.data)} directed edges, m range {M.data.min():.3g}-{M.data.max():.3g}", flush=True)

# ---------- figures: full and difference graphs ----------
projection = ccrs.EquidistantConic(central_longitude=coord[:, 0].mean(),
                                   central_latitude=coord[:, 1].mean())
fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300, subplot_kw={"projection": projection})
for ax, mode in zip(axs, ["Full", "Difference"]):
    v = Vis(ax, sp_digraph, projection=projection, edge_width=1, edge_alpha=1, edge_zorder=100,
            sample_pt_size=10, obs_node_size=5, sample_pt_color="black", cbar_font_size=6)
    v.draw_map()
    v.draw_samples()
    v.draw_edges(mode=mode)
    ax.set_title(f"{mode} graph", fontsize=9)
for ext in ("jpg", "pdf"):
    plt.savefig(f"{prefix}_frame_lambda{tag}.{ext}", bbox_inches="tight")
print(f"Saved: {prefix}_frame_lambda{tag}.jpg / .pdf / _edges.tsv / _nodes.tsv")
