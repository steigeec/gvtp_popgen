#!/usr/bin/env Rscript
# Initialize and run GENESPACE across the 6 species staged in genespace/peptide
# and genespace/bed, using the complete all-vs-all OrthoFinder run at
# genespace/orthofinder (a symlink: run_genespace() always reads OrthoFinder
# results from wd/orthofinder, whatever rawOrthofinderDir says).
#
# ploidy = 1 for every genome: none of the annotations are subgenome-phased.
# genespace_patches.R is sourced first (fix for run_mcscanx() crashing on
# synteny-poor genome pairs).
#
# Usage: Rscript run_genespace.R
# Conda env: genespace
#   built by annotation/setup_genespace.sh: conda env with r-base=4.3,
#   MCScanX and GENESPACE's Bioconductor/CRAN dependencies, then
#   devtools::install_github("jtlovell/GENESPACE") (GENESPACE 1.3.1)

BASE_DIR <- "/private/groups/meyerlab/esteiger"

library(GENESPACE)
source(file.path(BASE_DIR, "gvtp_popgen/inversion/genespace_patches.R"))

wd <- file.path(BASE_DIR, "inversions/synteny/genespace")

genomeIDs <- c(
  "Deinandra",
  "Helianthus_annuus",
  "Helianthus_anomalus",
  "Lactuca_sativa",
  "Arabidopsis_thaliana",
  "Solanum_lycopersicum"
)

gpar <- init_genespace(
  wd = wd,
  genomeIDs = genomeIDs,
  ploidy = rep(1, length(genomeIDs)),
  path2mcscanx = "/private/home/esteiger/miniconda3/envs/genespace/bin",
  rawOrthofinderDir = file.path(wd, "orthofinder"),
  nCores = 1
)

out <- run_genespace(gpar)

saveRDS(out, file.path(wd, "gsParam_final.rds"))
cat("Done. GENESPACE output under:", wd, "\n")
cat("  results/    - combined bed, orthogroups, blast-derived synteny hits\n")
cat("  riparian/   - genome-wide synteny plots (if plotting succeeded)\n")
cat("  dotplots/   - pairwise synteny dotplots\n")
