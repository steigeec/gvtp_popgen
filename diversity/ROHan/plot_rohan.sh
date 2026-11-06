#!/bin/bash
# =============================================================================
# plot_rohan.sh — Run rohan_plot.R after collect_rohan.sh completes
# =============================================================================

#SBATCH --job-name=ROHan_plot
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/diversity/ROH/ROHan/rohan_plot_%j.log

set -euo pipefail

Rscript /private/groups/meyerlab/esteiger/gvtp_popgen/diversity/ROHan/rohan_plot.R
