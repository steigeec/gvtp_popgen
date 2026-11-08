## Mantel test for isolation by distance at the individual level; companion to
## mantel_IBD.R (population centroids).
##
## Uses log(geographic distance), as in the population-level test (Rousset 1997).
## Samples sharing a coordinate (e.g. 15/16 CDFW samples, georeferenced at the
## site level) give zero distances, so they are collapsed to the highest-coverage
## individual first.
## Strauss (277 deduplicated individuals) is far denser than other clusters
## (1-15 each), so an all-individual test is dominated by within-Strauss pairs
## and between-cluster structure. Three tests are reported:
##   1. All deduplicated individuals, all clusters (descriptive)
##   2. Within Strauss only (fine-scale continuous IBD)
##   3. All clusters, with Strauss repeatedly subsampled to the size of the
##      next-largest cluster (sensitivity check)

BASE_DIR <- "/private/groups/meyerlab/esteiger"

.libPaths(file.path(BASE_DIR, "software/R_libs/"))
library(tidyverse)
library(geosphere)
library(vegan)

# === Load data ===
fst <- read.table("pairwise_Fst.txt", header = TRUE, sep = "\t",
                  colClasses = c("character", "character", "numeric", "numeric"))

meta <- read.table(file.path(BASE_DIR, "field_data.txt"),
                   header = TRUE, sep = "\t",
                   colClasses = c(SampleID = "character"))

bad_coords <- c("4548", "4549", "4550", "4551", "4552", "4553",
                "4554", "4555", "4556", "4557", "4559", "4561", "4562")

samples_used <- sort(setdiff(unique(c(fst$sample1, fst$sample2)), bad_coords))
meta_used <- meta %>% filter(SampleID %in% samples_used)

# === Collapse samples sharing a coordinate; keep highest-coverage individual ===
meta_used <- meta_used %>% mutate(coord = paste(Lat, Long))
keep <- meta_used %>%
  group_by(coord) %>%
  slice_max(Cov, n = 1, with_ties = FALSE) %>%
  ungroup()

cat(sprintf("Individuals before coordinate dedup: %d\n", nrow(meta_used)))
cat(sprintf("Individuals after coordinate dedup:  %d\n", nrow(keep)))
cat("\nCluster sizes after dedup:\n")
print(table(keep$Clusters))

fst <- fst %>%
  filter(sample1 %in% keep$SampleID, sample2 %in% keep$SampleID) %>%
  left_join(meta %>% select(SampleID, Lat1 = Lat, Long1 = Long, Cluster1 = Clusters), by = c("sample1" = "SampleID")) %>%
  left_join(meta %>% select(SampleID, Lat2 = Lat, Long2 = Long, Cluster2 = Clusters), by = c("sample2" = "SampleID")) %>%
  mutate(geo_dist_km = distHaversine(cbind(Long1, Lat1), cbind(Long2, Lat2)) / 1000,
         Fst_lin      = Fst_weighted / (1 - Fst_weighted),
         geo_log      = log(geo_dist_km))

stopifnot(all(is.finite(fst$geo_log)))  # confirms dedup removed all zero-distance pairs

n_all <- length(unique(c(fst$sample1, fst$sample2)))
cat(sprintf("\nDeduped Mantel dataset: %d individuals, %d pairs (expected %d)\n",
            n_all, nrow(fst), n_all * (n_all - 1) / 2))
stopifnot(nrow(fst) == n_all * (n_all - 1) / 2)

# === Helper: build dist objects + run Mantel from a pair-level data frame ===
run_mantel <- function(dat, permutations = 9999, seed = 42) {
  ids <- sort(unique(c(dat$sample1, dat$sample2)))
  n <- length(ids)
  fst_mat <- matrix(0, n, n, dimnames = list(ids, ids))
  geo_mat <- matrix(0, n, n, dimnames = list(ids, ids))
  for (i in seq_len(nrow(dat))) {
    a <- dat$sample1[i]; b <- dat$sample2[i]
    fst_mat[a, b] <- fst_mat[b, a] <- dat$Fst_lin[i]
    geo_mat[a, b] <- geo_mat[b, a] <- dat$geo_log[i]
  }
  set.seed(seed)
  res <- mantel(as.dist(fst_mat), as.dist(geo_mat), method = "pearson", permutations = permutations)
  list(n = n, n_pairs = nrow(dat), r = res$statistic, p = res$signif, perms = res$permutations)
}

# ============================================================
# Test 1: all deduped individuals, all clusters
# ============================================================
res_all <- run_mantel(fst)
cat("\n--- Test 1: all clusters, all deduped individuals ---\n")
cat(sprintf("n = %d individuals, %d pairs | Mantel r = %.4f | p = %s | permutations = %d\n",
            res_all$n, res_all$n_pairs, res_all$r,
            ifelse(res_all$p == 0, sprintf("< %.5g", 1 / (res_all$perms + 1)), signif(res_all$p, 3)),
            res_all$perms))

# ============================================================
# Test 2: within-Strauss only
# ============================================================
fst_strauss <- fst %>% filter(Cluster1 == "Strauss", Cluster2 == "Strauss")
res_strauss <- run_mantel(fst_strauss)
cat("\n--- Test 2: within-Strauss only ---\n")
cat(sprintf("n = %d individuals, %d pairs | Mantel r = %.4f | p = %s | permutations = %d\n",
            res_strauss$n, res_strauss$n_pairs, res_strauss$r,
            ifelse(res_strauss$p == 0, sprintf("< %.5g", 1 / (res_strauss$perms + 1)), signif(res_strauss$p, 3)),
            res_strauss$perms))

# ============================================================
# Test 3: all clusters, Strauss subsampled to match next-largest cluster,
#         repeated many times
# ============================================================
cluster_sizes <- table(keep$Clusters)
target_n <- max(cluster_sizes[names(cluster_sizes) != "Strauss"])
cat(sprintf("\n--- Test 3: Strauss subsampled to n = %d (size of largest non-Strauss cluster), %s ---\n",
            target_n, "100 iterations"))

strauss_ids  <- keep$SampleID[keep$Clusters == "Strauss"]
other_ids    <- keep$SampleID[keep$Clusters != "Strauss"]

n_iter <- 100
sub_r <- numeric(n_iter)
sub_p <- numeric(n_iter)
for (i in seq_len(n_iter)) {
  set.seed(1000 + i)
  strauss_sub <- sample(strauss_ids, target_n)
  keep_ids <- c(other_ids, strauss_sub)
  dat_sub <- fst %>% filter(sample1 %in% keep_ids, sample2 %in% keep_ids)
  r <- run_mantel(dat_sub, permutations = 999, seed = i)
  sub_r[i] <- r$r
  sub_p[i] <- r$p
}
cat(sprintf("Mantel r across %d subsamples: median = %.3f, range = [%.3f, %.3f]\n",
            n_iter, median(sub_r), min(sub_r), max(sub_r)))
cat(sprintf("Proportion of subsamples with p < 0.05: %.2f\n", mean(sub_p < 0.05)))

# ============================================================
# Plots
# ============================================================

make_mantel_plot <- function(dat, res, file, title_n) {
  fit <- lm(dat$Fst_lin ~ dat$geo_log)
  p_label <- if (res$p == 0) {
    sprintf("p < %.5g", 1 / (res$perms + 1))
  } else {
    paste0("p = ", signif(res$p, 3))
  }
  jpeg(file, width = 6, height = 5, units = "in", res = 300)
  plot(
    dat$geo_log, dat$Fst_lin,
    xlab = "ln(Geographic distance, km)",
    ylab = expression(Linearized~F[ST]),
    pch = 21, bg = adjustcolor("#009e73", alpha.f = 0.15), col = adjustcolor("black", alpha.f = 0.2),
    cex = 0.6
  )
  abline(fit, col = "black", lwd = 2)
  legend(
    "topleft", bty = "n",
    legend = c(
      paste0("Mantel r = ", round(res$r, 3)),
      p_label,
      title_n
    )
  )
  dev.off()
  cat("Saved", file, "\n")
}

make_mantel_plot(fst, res_all, "IBD_mantel_plot_individuals_all.jpg",
                  paste0("n = ", res_all$n, " individuals, ", res_all$n_pairs, " pairs"))
make_mantel_plot(fst_strauss, res_strauss, "IBD_mantel_plot_individuals_strauss.jpg",
                  paste0("n = ", res_strauss$n, " individuals, ", res_strauss$n_pairs, " pairs"))

jpeg("IBD_mantel_subsample_distribution.jpg", width = 5, height = 5, units = "in", res = 300)
hist(sub_r, breaks = 20, col = "#009e73", border = "white",
     main = "Mantel r across balanced Strauss subsamples",
     xlab = "Mantel r")
abline(v = res_all$r, col = "black", lwd = 2, lty = 2)
legend("topright", bty = "n", lty = 2, lwd = 2,
       legend = paste0("all-individuals r = ", round(res_all$r, 3)))
dev.off()
cat("Saved IBD_mantel_subsample_distribution.jpg\n")
