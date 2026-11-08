#!/usr/bin/env Rscript
# Convert ngsParalog likelihood ratios to p-values, plot diagnostics, and write
# postParalog.sites (sites passing the Bonferroni threshold plus sites that
# passed in an earlier run) and postParalog_MAF.sites (MAF > 0.05 subset).

BASE_DIR <- "/private/groups/meyerlab/esteiger"

library(tidyverse)
library(data.table)

PARA_DIR  <- file.path(BASE_DIR, "variant_calling/ngsParalog")
PREFIX    <- "noMAF_preParalog"
OUTDIR    <- file.path(BASE_DIR, "variant_calling/angsd")

# ---------------------------------------------------------------------------
# 1. Read ngsParalog output and compute p-values
# ---------------------------------------------------------------------------
cat("Reading ngsParalog LR file...\n")
lr <- fread(
  cmd = sprintf("awk '{print $1\"\t\"$2\"\t\"$5}' %s", file.path(PARA_DIR, paste0(PREFIX, ".lr"))),
  col.names    = c("chrom", "bp", "mlr"),
  colClasses   = list(factor = 1L)
) |> as.data.frame() |> unique()
lr$mlr <- as.numeric(as.character(lr$mlr))

lr$pval     <- 0.5 * pchisq(lr$mlr, df = 1, lower.tail = FALSE)
lr$pval.adj <- p.adjust(lr$pval, method = "bonferroni")

cat(sprintf("Total sites in LR file: %d\n", nrow(lr)))

# ---------------------------------------------------------------------------
# 2. Diagnostic plots
# ---------------------------------------------------------------------------
pdf(file.path(OUTDIR, paste0(PREFIX, "_paralog_diagnostics.pdf")), width = 10, height = 5)

# Histogram of raw p-values (uniform under null, excess near 0 for paralogs)
hist(lr$pval, breaks = 100, freq = FALSE,
     main = "Distribution of raw ngsParalog p-values",
     xlab = "p-value", col = "steelblue", border = "white")
abline(h = 1, col = "red", lty = 2)

# QQ-plot: observed vs expected -log10(p-values)
obs  <- sort(-log10(lr$pval), decreasing = TRUE)
expe <- -log10(ppoints(length(obs)))
plot(expe, obs, pch = 20, cex = 0.3,
     xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
     main = "QQ-plot: ngsParalog p-values")
abline(0, 1, col = "red")

dev.off()
cat(sprintf("Diagnostic plots saved to: %s\n",
            file.path(OUTDIR, paste0(PREFIX, "_paralog_diagnostics.pdf"))))

# ---------------------------------------------------------------------------
# 3. Retention summary at various thresholds
# ---------------------------------------------------------------------------
cat("\nngsParalog retention at various Bonferroni-adjusted thresholds:\n")
for (t in c(0.001, 0.005, 0.01, 0.05, 0.1)) {
  n <- sum(lr$pval.adj >= t, na.rm = TRUE)
  cat(sprintf("  pval.adj >= %-5s : %d sites retained (%.1f%%)\n",
              t, n, 100 * n / nrow(lr)))
}

#  pval.adj >= 0.001 : 37511492 sites retained (78.9%)
#  pval.adj >= 0.005 : 37225406 sites retained (78.3%)
#  pval.adj >= 0.01  : 37094848 sites retained (78.1%)
#  pval.adj >= 0.05  : 36772361 sites retained (77.4%)
#  pval.adj >= 0.1   : 36625486 sites retained (77.1%)
# I chose 0.05: few sites are lost at any of these thresholds, and no HWE
# p-value filter is applied first

# ---------------------------------------------------------------------------
# 4. Apply chosen threshold and export
# ---------------------------------------------------------------------------
# Set threshold here after inspecting the diagnostic plots above
PARALOG_THRESHOLD <- 0.05   # Bonferroni-adjusted p-value cutoff

cat(sprintf("\nApplying threshold: pval.adj >= %s\n", PARALOG_THRESHOLD))
sites_pass <- lr[lr$pval.adj >= PARALOG_THRESHOLD, c("chrom", "bp")]
cat(sprintf("  Sites retained from new LR: %d (%.1f%%)\n", nrow(sites_pass), 100 * nrow(sites_pass) / nrow(lr)))

# add sites that passed ngsParalog in an earlier run
prev_filtered <- fread(
  file.path(PARA_DIR, "previously-paralog-filtered.sites"),
  col.names   = c("chrom", "bp"),
  colClasses  = c("character", "integer")
) |> as.data.frame()
cat(sprintf("  Previously-filtered sites to add: %d\n", nrow(prev_filtered)))

sites_pass <- bind_rows(
  sites_pass |> mutate(chrom = as.character(chrom), bp = as.integer(bp)),
  prev_filtered
) |>
  unique() |>
  arrange(chrom, as.integer(bp))
cat(sprintf("  Total sites after merging: %d\n", nrow(sites_pass)))

outfile <- file.path(OUTDIR, "postParalog.sites")
write.table(sites_pass, outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("\nFiltered site list written to: %s\n", outfile))

# ---------------------------------------------------------------------------
# 5. MAF > 0.05 subset
# ---------------------------------------------------------------------------
maf <- fread(file.path(PARA_DIR, "preParalog_0.05maf.sites")) |> as.data.frame() |> unique()
names(maf)[1:2] <- c("chrom", "bp")

maf_paralog <- merge(maf, sites_pass, by = c("chrom", "bp")) |>
  arrange(chrom, as.integer(bp))

outfile <- file.path(OUTDIR, "postParalog_MAF.sites")
write.table(maf_paralog, outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Next: filterParalogMAF_prepNGSld.sh (LD pruning inputs)
