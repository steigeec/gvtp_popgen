#!/usr/bin/env Rscript
# Diagnostic: genomic inflation (lambda) and a small QQ plot from a raw GEMMA
# assoc.txt, without genomic-control rescaling.
#
# Usage: Rscript check_qq_lambda.R <assoc_txt_path> <out_png_path>

BASE_DIR <- "/private/groups/meyerlab/esteiger"

args <- commandArgs(trailingOnly = TRUE)
assoc_path <- args[1]
out_png    <- args[2]

.libPaths(file.path(BASE_DIR, "software/R_libs"))
library(data.table)
library(qqman)

gemma <- fread(assoc_path, header = TRUE)
print(paste("No rows in input:", nrow(gemma)))

chisq_obs <- qchisq(gemma$p_lrt, df = 1, lower.tail = FALSE)
lambda <- median(chisq_obs, na.rm = TRUE) / qchisq(0.5, df = 1)
print(paste("Lambda:", round(lambda, 3)))

dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
png(out_png, width = 1.5, height = 1.5, units = "in", res = 600, pointsize = 6)
par(mar = c(2.4, 2.4, 1.3, 0.3), mgp = c(1.3, 0.35, 0), tcl = -0.25,
    cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.75)
qq(gemma$p_lrt, main = paste0("lambda = ", round(lambda, 3)),
   xlab = "Expected -log10(p)", ylab = "Observed -log10(p)", pch = 20, cex = 0.4)
dev.off()
print(paste("Wrote", out_png))
