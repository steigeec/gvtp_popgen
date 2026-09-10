## AEM + RDA, HYBRID wind climatology variant: identical to wind_aem_rda.R
## except the directed wind graph comes from wind_support_asymmetric_hybrid.txt
## (WindNinja direction + raw HRRR speed; see build_hybrid_windrose.py and
## build_wind_support_hybrid.R) instead of wind_support_asymmetric.txt (pure
## HRRR direction+speed). Run after wind_aem_rda.R so both results exist for
## comparison -- see wind_aem_rda_results.txt vs wind_aem_rda_results_hybrid.txt.
##
## Rationale for the swap (not for re-deriving the allele-frequency pipeline,
## which is unchanged -- see wind_aem_rda.R's header for that): PTGC1
## validation showed WindNinja's terrain downscaling improves wind DIRECTION
## at an exposed cape station (median error 16 deg -> 10 deg) but does not
## recover HRRR's own speed underestimate there (28% -> 32% low, i.e. no
## better) -- so direction comes from WindNinja, speed stays from HRRR.

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(data.table)
library(vegan)
library(adespatial)
library(geosphere)

WIND_DIR <- "/private/groups/meyerlab/esteiger/wind"
ALLELE_FREQ_DIR <- "/private/groups/meyerlab/esteiger/wind/allele_freq"

maf_files <- c(
  CDFW             = "CDFW_maf.mafs.gz",
  Dangermond       = "Dangermond_maf.mafs.gz",
  Point_Conception = "Point_Conception_maf.mafs.gz",
  Strauss          = "Strauss_maf.mafs.gz",
  VAFB_Arguello    = "VAFB_Arguello_maf.mafs.gz",
  VAFB_Lion_Head   = "VAFB_Lion_Head_maf.mafs.gz"
)

t0 <- Sys.time()
pop_tables <- list()
for (pop in names(maf_files)) {
  dt <- fread(file.path(ALLELE_FREQ_DIR, maf_files[pop]),
              select = c("chromo", "position", "knownEM"))
  setnames(dt, "knownEM", pop)
  setkey(dt, chromo, position)
  pop_tables[[pop]] <- dt
  cat(sprintf("  %-17s %d sites [%.0fs]\n", pop, nrow(dt),
              as.numeric(Sys.time() - t0, units = "secs")))
}

wide <- Reduce(function(a, b) merge(a, b, by = c("chromo", "position"), all = TRUE), pop_tables)
for (pop in names(maf_files)) {
  set(wide, which(is.na(wide[[pop]])), pop, 0)
}
cat(sprintf("Merged wide table: %d sites x %d populations [%.0fs total]\n",
            nrow(wide), length(maf_files), as.numeric(Sys.time() - t0, units = "secs")))

mat <- as.matrix(wide[, names(maf_files), with = FALSE])
rownames(mat) <- paste(wide$chromo, wide$position, sep = ":")
Y_all <- t(mat)
storage.mode(Y_all) <- "double"

name_map <- c(
  "Point_Conception" = "conception", "CDFW" = "cdfw",
  "VAFB_Arguello"    = "arguello",   "Strauss" = "strauss",
  "VAFB_Lion_Head"   = "lionhead",   "Dangermond" = "Dangermond"
)
rownames(Y_all) <- name_map[rownames(Y_all)]

cat(sprintf("Allele frequency matrix: %d populations x %d loci\n", nrow(Y_all), ncol(Y_all)))

MIN_MAX_MAF <- 0.05

filter_sites <- function(Y, label) {
  max_freq <- apply(Y, 2, max)
  keep_maf <- max_freq > MIN_MAX_MAF
  cat(sprintf("  %s: %d / %d loci pass max-population-MAF > %.2f\n",
              label, sum(keep_maf), ncol(Y), MIN_MAX_MAF))
  Y <- Y[, keep_maf, drop = FALSE]
  keep_poly <- apply(Y, 2, function(x) var(x) > 0)
  cat(sprintf("  %s: %d / %d of those also polymorphic across these %d populations (kept)\n",
              label, sum(keep_poly), ncol(Y), nrow(Y)))
  Y[, keep_poly, drop = FALSE]
}

cat("Filtering loci (all 6 populations):\n")
Y_all <- filter_sites(Y_all, "6-pop")

Y_no_dang <- Y_all[rownames(Y_all) != "Dangermond", ]
cat("Filtering loci (5 populations, Dangermond excluded):\n")
Y_no_dang <- filter_sites(Y_no_dang, "5-pop")

field <- read.delim("/private/groups/meyerlab/esteiger/field_data.txt", header = TRUE, sep = "\t")
centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field, FUN = mean)
centroids$pop <- name_map[centroids$Clusters]
rownames(centroids) <- centroids$pop

coords_all <- as.matrix(centroids[rownames(Y_all), c("Long", "Lat")])
geo_all <- as.matrix(distm(coords_all, fun = distHaversine)) / 1000
dimnames(geo_all) <- list(rownames(Y_all), rownames(Y_all))

# HYBRID directed wind graph (WindNinja direction + HRRR speed) --
# see build_wind_support_hybrid.R.
wind_asym_full <- as.matrix(read.delim(file.path(WIND_DIR, "wind_support_asymmetric_hybrid.txt"), row.names = 1))

build_aem <- function(pops) {
  w <- wind_asym_full[pops, pops]
  n <- length(pops)
  edges <- expand.grid(from = pops, to = pops, stringsAsFactors = FALSE)
  edges <- edges[edges$from != edges$to, ]
  edges$weight <- mapply(function(i, j) w[i, j], edges$from, edges$to)
  se_mat <- matrix(0, nrow = n, ncol = nrow(edges), dimnames = list(pops, NULL))
  for (e in seq_len(nrow(edges))) {
    se_mat[edges$from[e], e] <- 1
    se_mat[edges$to[e], e] <- -1
  }
  edge_w <- edges$weight / max(edges$weight)
  ae <- aem(binary.mat = se_mat, weight = edge_w)
  list(vectors = ae$vectors, values = ae$values, edges = edges, se_mat = se_mat)
}

cat("\n=== Building AEM eigenvectors (6 populations, HYBRID wind) ===\n")
aem_all <- build_aem(rownames(Y_all))
cat("AEM eigenvalues (6 pop):\n"); print(round(aem_all$values, 4))

cat("\n=== Building AEM eigenvectors (5 populations, no Dangermond, HYBRID wind) ===\n")
aem_nd <- build_aem(rownames(Y_no_dang))
cat("AEM eigenvalues (5 pop):\n"); print(round(aem_nd$values, 4))

run_rda <- function(Y, geo, aem_vectors, label) {
  geo_pcoa <- cmdscale(as.dist(geo), k = 1)
  aem1 <- aem_vectors[, 1, drop = FALSE]
  colnames(aem1) <- "AEM1"
  df <- data.frame(geoPCoA1 = geo_pcoa[, 1], AEM1 = aem1[, 1])
  rownames(df) <- rownames(Y)

  mod <- rda(Y ~ AEM1 + Condition(geoPCoA1), data = df)
  set.seed(42)
  t0 <- Sys.time()
  test <- anova.cca(mod, permutations = 4999, by = NULL)
  cat(sprintf("[RDA + permutation test: %.0fs]\n", as.numeric(Sys.time() - t0, units = "secs")))
  cat(sprintf("\n--- RDA: %s ---\n", label))
  print(mod)
  cat("\nPermutation test (AEM1 | geoPCoA1):\n")
  print(test)
  list(model = mod, test = test)
}

cat("\n############ 6-population model (PRIMARY, HYBRID wind) ############\n")
res_all <- run_rda(Y_all, geo_all, aem_all$vectors, "6 populations, hybrid wind")

cat("\n############ 5-population model (Dangermond excluded, HYBRID wind) ############\n")
geo_nd <- geo_all[rownames(Y_no_dang), rownames(Y_no_dang)]
res_nd <- run_rda(Y_no_dang, geo_nd, aem_nd$vectors, "5 populations, no Dangermond, hybrid wind")

summary_txt <- file.path(WIND_DIR, "wind_aem_rda_results_hybrid.txt")
sink(summary_txt)
cat("AEM + RDA directional wind connectivity test -- HYBRID wind (WindNinja direction + HRRR speed)\n")
cat("=================================================================================\n\n")
cat("PRIMARY result: 6 populations\n")
print(res_all$test)
cat("\n\nRobustness check: 5 populations, Dangermond excluded\n")
print(res_nd$test)
sink()
cat("\nSaved summary ->", summary_txt, "\n")

# ---------- 6. RDA triplot: visualize the primary 6-population, hybrid-wind
#              result (site scores on the single constrained "wind" axis vs.
#              the strongest residual axis, with the AEM1/wind constraint
#              shown as a biplot arrow -- geoPCoA1 is a Condition() term, so
#              it is partialled out of the ordination and does not appear). ----------

saveRDS(res_all$model, file.path(WIND_DIR, "rda_model_hybrid_6pop.rds"))
saveRDS(res_nd$model, file.path(WIND_DIR, "rda_model_hybrid_5pop.rds"))

pop_labels <- c(
  lionhead = "Lion's Head", arguello = "Pt. Arguello", strauss = "Tranquillon",
  Dangermond = "Dangermond", conception = "Pt. Conception", cdfw = "Gaviota"
)

mod <- res_all$model
xy <- scores(mod, display = "sites", scaling = 2, choices = c(1, 2))
cat("\nOrdination axes used for the triplot:", paste(colnames(xy), collapse = ", "), "\n")

jpeg(file.path(WIND_DIR, "wind_aem_rda_triplot.jpg"), width = 6, height = 5.5, units = "in", res = 300)
plot(xy[, 1], xy[, 2], type = "n",
     xlab = sprintf("RDA1 (wind, AEM1 | geography) -- %.1f%% of total variance",
                     100 * mod$CCA$eig[1] / mod$tot.chi),
     ylab = sprintf("PC1 (residual) -- %.1f%% of total variance",
                     100 * mod$CA$eig[1] / mod$tot.chi),
     main = "RDA: allele frequencies ~ wind (AEM1) | geography\n(6 populations, hybrid wind climatology)")
abline(h = 0, v = 0, col = "grey80", lty = 2)
points(xy[, 1], xy[, 2], pch = 19, col = "#0072b2", cex = 1.6)
text(xy[, 1], xy[, 2], labels = pop_labels[rownames(xy)], pos = 3, cex = 0.9, font = 2)

bp <- scores(mod, display = "bp", scaling = 2)
bp1 <- bp[1, c(1, 2)]  # the single constraining variable, AEM1
arrow_scale <- 0.8 * max(abs(xy)) / max(abs(bp1))
arrows(0, 0, bp1[1] * arrow_scale, bp1[2] * arrow_scale, length = 0.1, col = "firebrick", lwd = 2)
text(bp1[1] * arrow_scale, bp1[2] * arrow_scale, labels = "wind (AEM1)", col = "firebrick", pos = 4, font = 2)

mtext(sprintf("Permutation test: F = %.2f, p = %.2f (n.s.)",
              res_all$test$F[1], res_all$test$"Pr(>F)"[1]), side = 3, line = 0.3, cex = 0.85)
dev.off()
cat("Saved ->", file.path(WIND_DIR, "wind_aem_rda_triplot.jpg"), "\n")
