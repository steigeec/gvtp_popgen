#!/usr/bin/env Rscript
# Compares observed lambda and Bonferroni hit count on a SNP subset against a
# phenotype-permutation null; reports empirical p-values and a histogram figure.
#
# Usage: Rscript summarize_permutation_null.R <var_name> <n_perm> <out_dir>

BASE_DIR <- "/private/groups/meyerlab/esteiger"

args     <- commandArgs(trailingOnly = TRUE)
var_name <- args[1]
n_perm   <- as.integer(args[2])
out_dir  <- args[3]

.libPaths(file.path(BASE_DIR, "software/R_libs"))
library(data.table)

get_pvec <- function(f) fread(f, header = TRUE, select = "p_lrt")$p_lrt

lambda_of <- function(p) median(qchisq(p, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, df = 1)
hits_of   <- function(p, thresh) sum(p < thresh, na.rm = TRUE)

obs_p <- get_pvec(paste0("output/", var_name, "_subset_observed.assoc.txt"))
n_snps <- length(obs_p)
thresh <- 0.05 / n_snps

lambda_obs <- lambda_of(obs_p)
hits_obs   <- hits_of(obs_p, thresh)

lambda_perm <- numeric(n_perm)
hits_perm   <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  p <- get_pvec(paste0("output/", var_name, "_perm_", i, ".assoc.txt"))
  lambda_perm[i] <- lambda_of(p)
  hits_perm[i]   <- hits_of(p, thresh)
}

p_lambda <- (sum(lambda_perm >= lambda_obs) + 1) / (n_perm + 1)
p_hits   <- (sum(hits_perm   >= hits_obs)   + 1) / (n_perm + 1)

cat("=== Permutation null summary for", var_name, "===\n")
cat("N SNPs in subset:", n_snps, "\n")
cat("Bonferroni threshold (subset-scaled):", thresh, "\n")
cat("Observed lambda:", round(lambda_obs, 3), "\n")
cat("Permutation lambda: mean =", round(mean(lambda_perm), 3),
    " range =", round(min(lambda_perm), 3), "-", round(max(lambda_perm), 3), "\n")
cat("Empirical p (lambda):", p_lambda, "\n")
cat("Observed hits:", hits_obs, "\n")
cat("Permutation hits: mean =", round(mean(hits_perm), 2),
    " range =", min(hits_perm), "-", max(hits_perm), "\n")
cat("Empirical p (hits):", p_hits, "\n")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
png(file.path(out_dir, paste0("permutation_null_", var_name, ".png")),
    width = 5, height = 2.5, units = "in", res = 300, pointsize = 9)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
hist(lambda_perm, main = "Lambda: permutation null", xlab = "lambda",
     xlim = range(c(lambda_perm, lambda_obs)))
abline(v = lambda_obs, col = "red", lwd = 2)
hist(hits_perm, main = "Hits: permutation null", xlab = "N Bonferroni hits",
     xlim = range(c(hits_perm, hits_obs)))
abline(v = hits_obs, col = "red", lwd = 2)
dev.off()
cat("Wrote", file.path(out_dir, paste0("permutation_null_", var_name, ".png")), "\n")
