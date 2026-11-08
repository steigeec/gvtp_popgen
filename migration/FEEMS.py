#!/usr/bin/env python
"""
Fit FEEMS for one sampling year.

Usage: python FEEMS.py <plink_prefix>
Conda env: feems_e
  conda create -n feems_e -c conda-forge -c bioconda python=3.8.10 \
    numpy=1.22.3 scipy=1.5.0 scikit-learn=0.23.1 matplotlib=3.2.2 \
    networkx=2.4 pyproj=2.6.1 shapely=1.7.1 fiona cartopy=0.18.0 \
    suitesparse=5.7.2 scikit-sparse=0.4.4 pandas-plink msprime=1.0.0 \
    statsmodels=0.14.0
  pip install git+https://github.com/NovembreLab/feems@c08a179   # FEEMS 1.0.1
"""

#Import required packages
import numpy as np
import pkg_resources
import scipy as sp
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
import statsmodels.api as sm
import csv
import sys
 
# viz
import matplotlib.pyplot as plt
from matplotlib import gridspec
import cartopy.crs as ccrs
from pyproj import Proj
from cartopy.crs import Projection
# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz
from feems.cross_validation import run_cv 
# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"

# If using sys.argv:
if len(sys.argv) < 2:
    print("Usage: python FEEMS.py <plink_prefix>")
    sys.exit(1)

filename = sys.argv[1]

# read the genotype data, where "HET" is the basename
(bim, fam, G) = read_plink(filename)
 
# Count total missing genotypes across all genomic positions and individuals
total_missing_genotypes_dask = np.sum(np.isnan(G))
total_missing_genotypes = total_missing_genotypes_dask.compute()
print("Total Missing Genotypes:", total_missing_genotypes)   # 0 
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)
 
# Filter genotypes matrix based on variable sites
variable_sites = np.var(genotypes, axis=0) > 0  # Sites with variance > 0 are variable
genotypes = genotypes[:, variable_sites]
print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))
# n_samples=277, n_snps=2,588,155
 
# setup graph
coord_file = f"{filename}.coord"
coord = np.loadtxt(coord_file)  
outer = np.loadtxt("outer.csv")  
grid = np.loadtxt("grid.csv")  
edges = np.loadtxt("edges.csv", dtype=int)  
 
# Set up the spatialgrph object
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)


#set up the projection for this dataset
#isea_proj = Proj("+proj=isea +lon_0=-120.4 +lat_0=34.65")
#projection = Projection(isea_proj)
projection = ccrs.EquidistantConic(central_longitude=coord[:, 0].mean(), central_latitude=coord[:, 1].mean())

##### CROSS VALIDATION ############

# define grids
# reverse the order of lambdas and alphas for warmstart
# geomspace: Return numbers spaced evenly on a log scale (a geometric progression).  Start number, stop number, and number of intervals. 
# I include larger numbers than in example because of error below, which suggests unstable optimization:
#sksparse.cholmod.CholmodNotPositiveDefiniteError: ... not positive definite (code 1)
lamb_grid = np.geomspace(1e-6, 1e4, 2000)[::-1]
 
# run cross-validation
cv_err = run_cv(sp_graph, lamb_grid, n_folds=sp_graph.n_observed_nodes, factr=1e10)
 
# average over folds
mean_cv_err = np.mean(cv_err, axis=0)
 
# argmin of cv error
lamb_cv = float(lamb_grid[np.argmin(mean_cv_err)])
 
# Plot CV error:
fig, ax = plt.subplots(dpi=300)
ax.plot(np.log10(lamb_grid), mean_cv_err, ".");
ax.set_xlabel("log10(lambda)");
ax.set_ylabel("L2 CV Error");
ax.axvline(np.log10(lamb_cv), color = "orange")
plt.savefig('CV-error.jpg')


# Lets take a look at a summary of CV:
# figure params
projection = ccrs.EquidistantConic(central_longitude=coord[:, 0].mean(), central_latitude=coord[:, 1].mean())
title_loc = "left"
title_pad = "-10"
title_fontsize = 12
edge_width = 1
edge_alpha = 1
edge_zorder = 3
obs_node_size = 3
obs_node_linewidth = .4
cbar_font_size = 8
cbar_ticklabelsize = 8
cbar_orientation = "horizontal"
 
# figure setup
fig = plt.figure(figsize=(8, 6), dpi=300)
spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig, wspace=0.0, hspace=0.0)
 
# axis 00
ax_00 = fig.add_subplot(spec[0, 0], projection=projection)
ax_00.set_title("A", loc=title_loc, pad=title_pad, fontdict={"fontsize": title_fontsize})
sp_graph.fit(float(lamb_grid[0]))
v = Viz(ax_00, sp_graph, projection=projection, edge_width=edge_width,
        edge_alpha=1, edge_zorder=100, sample_pt_size=20,
        obs_node_size=obs_node_size, sample_pt_color="black",
        cbar_font_size=10)
v.draw_map()
v.draw_samples() 
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
ax_00.text(.2, .85, "lambda={:.5f}\ncv l2 error={:.5f}".format(lamb_grid[0], mean_cv_err[0, 0]),
           fontdict={"fontsize": 4}, transform = ax_00.transAxes)
ax.set_aspect('equal')
 
# axis 10
ax_10 = fig.add_subplot(spec[1, 0], projection=projection)
ax_10.set_title("B", loc=title_loc, pad=title_pad, fontdict={"fontsize": title_fontsize})
sp_graph.fit(float(lamb_grid[3]))
v = Viz(ax_10, sp_graph, projection=projection, edge_width=edge_width,
        edge_alpha=1, edge_zorder=100, sample_pt_size=20,
        obs_node_size=obs_node_size, sample_pt_color="black",
        cbar_font_size=10)
v.draw_map()
v.draw_samples() 
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
ax_10.text(.2, .85, "lambda={:.5f}\ncv l2 error={:.5f}".format(lamb_grid[3], mean_cv_err[3, 0]),
           fontdict={"fontsize": 4}, transform = ax_10.transAxes)
ax.set_aspect('equal')
 
# axis 01
ax_01 = fig.add_subplot(spec[0, 1], projection=projection)
ax_01.set_title("C", loc=title_loc, pad=title_pad, fontdict={"fontsize": title_fontsize})
sp_graph.fit(float(lamb_cv))
v = Viz(ax_01, sp_graph, projection=projection, edge_width=edge_width,
        edge_alpha=1, edge_zorder=100, sample_pt_size=20,
        obs_node_size=obs_node_size, sample_pt_color="black",
        cbar_font_size=10)
v.draw_map()
v.draw_samples() 
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
ax_01.text(.2, .85, "lambda={:.5f}\ncv l2 error={:.5f}".format(lamb_cv, mean_cv_err[np.argmin(mean_cv_err), 0]),
           fontdict={"fontsize": 4}, transform = ax_01.transAxes)
ax.set_aspect('equal')
 
# axis 11
ax_11 = fig.add_subplot(spec[1, 1], projection=projection)
ax_11.set_title("D", loc=title_loc, pad=title_pad, fontdict={"fontsize": title_fontsize})
sp_graph.fit(float(lamb_grid[10]))
v = Viz(ax_11, sp_graph, projection=projection, edge_width=edge_width,
        edge_alpha=1, edge_zorder=100, sample_pt_size=20,
        obs_node_size=obs_node_size, sample_pt_color="black",
        cbar_font_size=10)
v.draw_map()
v.draw_samples() 
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
v.cbar_font_size = cbar_font_size
v.cbar_orientation = cbar_orientation
v.cbar_ticklabelsize = cbar_ticklabelsize
v.draw_edge_colorbar()
ax_11.text(.2, .85, "lambda={:.5f}\ncv l2 error={:.5f}".format(lamb_grid[10], mean_cv_err[10, 0]),
           fontdict={"fontsize": 4}, transform = ax_11.transAxes)
ax.set_aspect('equal')
plt.savefig('CV-summary.jpg', bbox_inches='tight')
 
# Lets zoom in on the fit that has the lowest error on the CV curve:
 
# re-fit
sp_graph.fit(lamb_cv)
 
# Create a figure with a specific size
fig = plt.figure(figsize=(8, 6), dpi=300)             # Adjust the size as needed
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=1,    #increased edge size relative to original
        edge_alpha=1, edge_zorder=100, sample_pt_size=20,
        obs_node_size=7.5, sample_pt_color="black",
        cbar_font_size=10)
v.draw_map()
v.draw_samples() 
v.draw_edges(use_weights=True)  
v.draw_obs_nodes(use_ids=False)
v.draw_edge_colorbar()
# Set aspect ratio to 'equal' to fix axes proportions
ax.set_aspect('equal')
# Save the figure without extra white margins
plt.savefig(f'{filename}_LowestCVerror.jpg', bbox_inches='tight')

