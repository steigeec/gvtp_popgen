## Mantel test for isolation by distance among Deinandra increscens populations:
## Fst/(1-Fst) against ln(geographic distance) (Rousset 1997).
## Dangermond's pairwise Fst is inflated by low coverage (three of four
## individuals at 0.9-2.9x), so three versions are tested:
##   primary    5 populations (Dangermond excluded), full-sample Fst
##   all6_full  all 6 populations, full-sample Fst
##   all6_cov   all 6 populations, coverage-matched n = 4 Fst (mean of 3 draws)
## Output (OUT_DIR): IBD_mantel_results.txt, IBD_mantel_plot.jpg

BASE_DIR <- "/private/groups/meyerlab/esteiger"
.libPaths(file.path(BASE_DIR, "software/R_libs/"))
library(vegan)
library(geosphere)

OUT_DIR <- file.path(BASE_DIR, "structure/Fst/populations/IBD")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

FULL_FST <- file.path(BASE_DIR, "structure/Fst/populations/pairwise_auto/pairwise_Fst_populations.txt")
COVMATCH_FST <- file.path(BASE_DIR, "wind/fst_downsampled_covmatch/pairwise_Fst_downsampled.txt")

# ---------- 1. Compute population centroids from field data ----------

# 13 individuals have untrustworthy coordinates (Location "Strauss" but
# clustered with Point_Conception or VAFB_Lion_Head). Excluded from centroids.
BAD_COORD_IDS <- c("4548", "4549", "4550", "4551", "4552", "4553", "4554",
                    "4555", "4556", "4557", "4559", "4561", "4562")
field <- read.delim(
  file.path(BASE_DIR, "field_data.txt"),
  header = TRUE, sep = "\t"
)
field <- field[!as.character(field$SampleID) %in% BAD_COORD_IDS, ]

centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field, FUN = mean)

# Cluster names to the population names in the Fst output
name_map <- c(
  "Point_Conception" = "conception",
  "CDFW"             = "cdfw",
  "VAFB_Arguello"    = "arguello",
  "Strauss"          = "strauss",
  "VAFB_Lion_Head"   = "lionhead",
  "Dangermond"       = "Dangermond"
)
centroids$pop <- name_map[centroids$Clusters]
rownames(centroids) <- centroids$pop

cat("Population centroids:\n")
print(centroids[, c("pop", "Lat", "Long")])

# ---------- 2. Read pairwise Fst ----------

# Symmetric weighted-Fst matrix. Replicate pairs in the coverage-matched file
# are averaged.
read_fst <- function(file) {
  fst_df <- read.delim(file, header = TRUE, sep = "\t")
  fst_df <- aggregate(Fst_weighted ~ pop1 + pop2, data = fst_df, FUN = mean)
  pops <- sort(unique(c(fst_df$pop1, fst_df$pop2)))
  n <- length(pops)
  if (nrow(fst_df) < n * (n - 1) / 2) warning("Fst matrix is incomplete -- some pairs are missing: ", file)
  fst_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(pops, pops))
  fst_mat[cbind(fst_df$pop1, fst_df$pop2)] <- fst_df$Fst_weighted
  fst_mat[cbind(fst_df$pop2, fst_df$pop1)] <- fst_df$Fst_weighted
  fst_mat
}

fst_full <- read_fst(FULL_FST)
fst_cov <- read_fst(COVMATCH_FST)
all_pops <- rownames(fst_full)
stopifnot(setequal(all_pops, rownames(fst_cov)))

# ---------- 3. Geographic distance (km) and linearization (Rousset 1997) ----------

coords <- centroids[all_pops, c("Long", "Lat")]
geo_mat <- as.matrix(distm(coords, fun = distHaversine)) / 1000
dimnames(geo_mat) <- list(all_pops, all_pops)

cat("\nGeographic distances (km):\n")
print(round(geo_mat, 2))

geo_lin <- log(geo_mat)
diag(geo_lin) <- 0

# ---------- 4. Mantel tests ----------

run_mantel <- function(fst_mat, pops) {
  fst_lin <- fst_mat[pops, pops] / (1 - fst_mat[pops, pops])
  set.seed(42)
  res <- mantel(as.dist(fst_lin), as.dist(geo_lin[pops, pops]), method = "pearson", permutations = 9999)
  list(r = res$statistic, p = res$signif, perms = res$permutations,
       n_pops = length(pops), n_pairs = length(pops) * (length(pops) - 1) / 2, fst_lin = fst_lin)
}

pops5 <- setdiff(all_pops, "Dangermond")
tests <- list(
  primary   = c(run_mantel(fst_full, pops5),    fst = "full-sample", pops = "5 (Dangermond excluded)"),
  all6_full = c(run_mantel(fst_full, all_pops), fst = "full-sample", pops = "6"),
  all6_cov  = c(run_mantel(fst_cov, all_pops),  fst = "coverage-matched n = 4", pops = "6")
)

# With 5-6 populations there are only 120-720 permutations, so vegan
# enumerates them all (exact test) instead of drawing 9999.
cat("\n--- Mantel test results (Fst/(1-Fst) ~ ln km; one-sided, exact permutation test) ---\n")
summary_df <- do.call(rbind, lapply(names(tests), function(k) {
  t <- tests[[k]]
  data.frame(test = k, populations = t$pops, fst = t$fst, n_pairs = t$n_pairs,
             mantel_r = round(t$r, 4), p = round(t$p, 4), permutations = t$perms)
}))
print(summary_df, row.names = FALSE)
write.table(summary_df, file.path(OUT_DIR, "IBD_mantel_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- 5. Plot IBD ----------
# All 15 pairs, full-sample Fst. Dangermond pairs are hollow and outlined in
# the partner population's colour. Fit line and statistics from the primary test.

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

# One row per pair. Fill is the northern population, border the southern.
idx <- which(lower.tri(geo_lin), arr.ind = TRUE)
pair_df <- data.frame(p1 = all_pops[idx[, 1]], p2 = all_pops[idx[, 2]])
pair_df$geo <- geo_lin[idx]
pair_df$fst <- (fst_full / (1 - fst_full))[idx]
north <- ifelse(centroids[pair_df$p1, "Lat"] >= centroids[pair_df$p2, "Lat"], pair_df$p1, pair_df$p2)
south <- ifelse(north == pair_df$p1, pair_df$p2, pair_df$p1)
pair_df$fill <- popcolors[pop_to_cluster[north]]
pair_df$border <- popcolors[pop_to_cluster[south]]
dang <- pair_df$p1 == "Dangermond" | pair_df$p2 == "Dangermond"
partner <- ifelse(pair_df$p1 == "Dangermond", pair_df$p2, pair_df$p1)
pair_df$border[dang] <- popcolors[pop_to_cluster[partner[dang]]]

p_label <- function(t) {
  if (t$p * (t$perms + 1) <= 1) sprintf("p < %.5g", 1 / (t$perms + 1)) else sprintf("p = %.3g", t$p)
}

plot_file <- file.path(OUT_DIR, "IBD_mantel_plot.jpg")
jpeg(plot_file, width = 5, height = 5, units = "in", res = 300)
plot(
  pair_df$geo, pair_df$fst,
  xlab = "ln(Geographic distance, km)",
  ylab = expression(F[ST] / (1 - F[ST])),
  type = "n"
)
points(pair_df$geo[!dang], pair_df$fst[!dang], pch = 21,
       bg = pair_df$fill[!dang], col = pair_df$border[!dang], cex = 1.5, lwd = 3)
points(pair_df$geo[dang], pair_df$fst[dang], pch = 21,
       bg = "white", col = pair_df$border[dang], cex = 1.5, lwd = 3)
fit <- lm(fst ~ geo, data = pair_df[!dang, ])
abline(fit, col = "black", lwd = 2)
legend(
  "topleft", bty = "n",
  legend = c(
    paste0("Mantel r = ", round(tests$primary$r, 3)),
    p_label(tests$primary),
    "(5 populations; open symbols:",
    " Dangermond pairs, not tested;",
    " outline = paired population)"
  )
)
dev.off()
cat("\nResults saved to", OUT_DIR, "(IBD_mantel_results.txt, IBD_mantel_plot.jpg)\n")
