## Standalone replot of the primary 6-population hybrid-wind RDA triplot,
## loading the already-fitted model (saved by wind_aem_rda_hybrid.R) instead
## of rerunning the ~1.5-2 hour fit + permutation test. Fixes layout issues
## in the first version: title colliding with the Dangermond point, the
## wind(AEM1) arrow line running through the Lion's Head label, and Pt.
## Arguello's label clipped at the plot edge.

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(vegan)

WIND_DIR <- "/private/groups/meyerlab/esteiger/wind"
mod <- readRDS(file.path(WIND_DIR, "rda_model_hybrid_6pop.rds"))

# F/p from the already-completed permutation test (wind_aem_rda_results_hybrid.txt) --
# deterministic given seed=42, not recomputed here.
F_stat <- 0.5512
p_val <- 0.7222

pop_labels <- c(
  lionhead = "Lion's Head", arguello = "Pt. Arguello", strauss = "Tranquillon",
  Dangermond = "Dangermond", conception = "Pt. Conception", cdfw = "Gaviota"
)

xy <- scores(mod, display = "sites", scaling = 2, choices = c(1, 2))
bp <- scores(mod, display = "bp", scaling = 2)
bp1 <- bp[1, c(1, 2)]

xr <- extendrange(xy[, 1], f = 0.28)
yr <- extendrange(xy[, 2], f = 0.28)
arrow_scale <- 0.55 * max(abs(xy)) / max(abs(bp1))

# Per-point label offsets, chosen by hand to avoid collisions given this
# specific point layout (Dangermond sits highest -> label below; Tranquillon
# and Pt. Arguello sit close together at bottom-right -> stagger left/right;
# Lion's Head sits near the wind-arrow's line -> label below, away from it).
label_pos <- c(lionhead = 1, arguello = 4, strauss = 2,
                Dangermond = 1, conception = 3, cdfw = 1)

jpeg(file.path(WIND_DIR, "wind_aem_rda_triplot.jpg"), width = 6.5, height = 6, units = "in", res = 300)
par(mar = c(4.5, 4.5, 4.5, 2))
plot(xy[, 1], xy[, 2], type = "n", xlim = xr, ylim = yr,
     xlab = sprintf("RDA1 (wind, AEM1 | geography) -- %.1f%% of total variance",
                     100 * mod$CCA$eig[1] / mod$tot.chi),
     ylab = sprintf("PC1 (residual) -- %.1f%% of total variance",
                     100 * mod$CA$eig[1] / mod$tot.chi),
     main = "")
title(main = "RDA: allele frequencies ~ wind (AEM1) | geography\n(6 populations, hybrid wind climatology)",
      line = 2.6, cex.main = 1.15)
mtext(sprintf("Permutation test: F = %.2f, p = %.2f (n.s.)", F_stat, p_val),
      side = 3, line = 0.5, cex = 0.9)
abline(h = 0, v = 0, col = "grey80", lty = 2)

arrows(0, 0, bp1[1] * arrow_scale, bp1[2] * arrow_scale, length = 0.1, col = "firebrick", lwd = 2)
text(bp1[1] * arrow_scale, bp1[2] * arrow_scale, labels = "wind (AEM1)",
     col = "firebrick", pos = 1, cex = 0.95, font = 2)

points(xy[, 1], xy[, 2], pch = 19, col = "#0072b2", cex = 1.6)
text(xy[, 1], xy[, 2], labels = pop_labels[rownames(xy)],
     pos = label_pos[rownames(xy)], cex = 0.95, font = 2)

dev.off()
cat("Saved ->", file.path(WIND_DIR, "wind_aem_rda_triplot.jpg"), "\n")
