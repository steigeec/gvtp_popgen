## Mantel test for isolation by distance (IBD)
## Tests whether pairwise Fst correlates with geographic distance
## among Deinandra increscens populations.
##
## Uses Rousset (1997) linearization: Fst/(1-Fst) vs log(geographic distance)
## for 2D stepping-stone expectations.
##
## Possible reasons Fst may NOT correlate with geographic distance:
##  - Recent colonization or founder effects: genetic similarity reflects
##    colonization history rather than current spatial arrangement.
##  - Barriers to gene flow: physical or habitat discontinuities (ridges,
##    highways, unsuitable habitat) can isolate nearby populations.
##  - Genetic drift in small populations: small Ne causes rapid divergence
##    regardless of distance, inflating Fst for those comparisons.
##  - Human-mediated movement: restoration plantings or accidental transport
##    can make distant populations unexpectedly similar.
##  - Divergent selection / local adaptation: different microhabitats
##    (coastal vs. inland, soil type) can drive excess differentiation.
##  - Non-equilibrium: IBD requires drift-migration balance; young or
##    recently bottlenecked populations may not have reached equilibrium.

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(vegan)
library(geosphere)

# ---------- 1. Compute population centroids from field data ----------

field <- read.delim(
  "/private/groups/meyerlab/esteiger/field_data.txt",
  header = TRUE, sep = "\t"
)

centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field, FUN = mean)

# Map cluster names to the population names used in the Fst output
name_map <- c(
  "Point_Conception" = "conception",
  "CDFW"             = "cdfw",
  "VAFB_Arguello"    = "arguello",
  "Strauss"          = "strauss",
  "VAFB_Lion_Head"   = "lionhead",
  "Dangermond"       = "Dangermond"
)
centroids$pop <- name_map[centroids$Clusters]

cat("Population centroids:\n")
print(centroids[, c("pop", "Lat", "Long")])

# ---------- 2. Read pairwise Fst ----------

fst_df <- read.delim(
  "/private/groups/meyerlab/esteiger/structure/Fst/populations/pairwise_auto/pairwise_Fst_populations.txt",
  header = TRUE, sep = "\t"
)

pops <- sort(unique(c(fst_df$pop1, fst_df$pop2)))
n <- length(pops)

cat("\nPopulations in Fst file:", paste(pops, collapse = ", "), "\n")
cat("Expected pairs:", n * (n - 1) / 2, " | Found:", nrow(fst_df), "\n")

if (nrow(fst_df) < n * (n - 1) / 2) {
  warning("Fst matrix is incomplete -- some pairs are missing.")
}

# Build symmetric Fst matrix (weighted Fst)
fst_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(pops, pops))
for (i in seq_len(nrow(fst_df))) {
  p1 <- fst_df$pop1[i]
  p2 <- fst_df$pop2[i]
  fst_mat[p1, p2] <- fst_df$Fst_weighted[i]
  fst_mat[p2, p1] <- fst_df$Fst_weighted[i]
}

# ---------- 3. Build geographic distance matrix (km) ----------

centroids <- centroids[match(pops, centroids$pop), ]
coords <- centroids[, c("Long", "Lat")]
rownames(coords) <- centroids$pop

geo_mat <- as.matrix(distm(coords, fun = distHaversine)) / 1000
dimnames(geo_mat) <- list(pops, pops)

cat("\nGeographic distances (km):\n")
print(round(geo_mat, 2))

# ---------- 4. Linearize (Rousset 1997) ----------

fst_lin <- fst_mat / (1 - fst_mat)
geo_lin <- log(geo_mat)
diag(geo_lin) <- 0

# ---------- 5. Mantel test ----------

fst_dist <- as.dist(fst_lin)
geo_dist <- as.dist(geo_lin)

set.seed(42)
mantel_result <- mantel(fst_dist, geo_dist, method = "pearson", permutations = 9999)

cat("\n--- Mantel test results ---\n")
cat("Mantel r:", round(mantel_result$statistic, 4), "\n")
cat("p-value: ", round(mantel_result$signif, 4), "\n")
cat("Permutations:", mantel_result$permutations, "\n")

# ---------- 6. Plot IBD ----------

popcolors <- c(
  "CDFW"           = "#e69f00",
  "Strauss"        = "#009e73",
  "Point_Conception" = "#0072b2",
  "VAFB_Lion_Head"  = "#f0e442",
  "VAFB_Arguello"   = "#d55e00",
  "Dangermond"      = "#56b4e9"
)

pop_to_cluster <- c(
  "cdfw"       = "CDFW",
  "strauss"    = "Strauss",
  "conception" = "Point_Conception",
  "lionhead"   = "VAFB_Lion_Head",
  "arguello"   = "VAFB_Arguello",
  "Dangermond" = "Dangermond"
)

pop_lat <- setNames(centroids$Lat, centroids$pop)

pair_fill <- character(0)
pair_border <- character(0)
for (i in seq_len(nrow(fst_df))) {
  p1 <- fst_df$pop1[i]
  p2 <- fst_df$pop2[i]
  if (pop_lat[p1] >= pop_lat[p2]) {
    north <- p1; south <- p2
  } else {
    north <- p2; south <- p1
  }
  pair_fill   <- c(pair_fill,   popcolors[pop_to_cluster[north]])
  pair_border <- c(pair_border, popcolors[pop_to_cluster[south]])
}

geo_vec <- as.vector(geo_dist)
fst_vec <- as.vector(fst_dist)

jpeg("IBD_mantel_plot.jpg", width = 5, height = 5, units = "in", res = 300)
plot(
  geo_vec, fst_vec,
  xlab = "ln(Geographic distance, km)",
  ylab = expression(F[ST] / (1 - F[ST])),
  pch = 21, bg = pair_fill, col = pair_border,
  cex = 1.5, lwd = 3
)
fit <- lm(fst_vec ~ geo_vec)
abline(fit, col = "black", lwd = 2)
legend(
  "topleft", bty = "n",
  legend = c(
    paste0("Mantel r = ", round(mantel_result$statistic, 3)),
    paste0("p = ", round(mantel_result$signif, 4))
  )
)
dev.off()
cat("\nPlot saved to IBD_mantel_plot.jpg\n")
