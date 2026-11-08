## Edge-by-edge test of whether FRAME gene-flow asymmetry matches windscape
## wind-flow asymmetry on the FEEMS grid (grid.csv, edges.csv).
## Gene-flow asymmetry is ln(rate a to b / rate b to a) from the geneflow_from/_to
## columns, i.e. forward gene flow, the reverse of FRAME's lineage direction.
## Wind asymmetry is ln(hours b to a / hours a to b) from least-cost travel times
## on the 100 m WindNinja field (trans = 1).
## Null: the wind field rotated by N_ROT evenly spaced angles; statistic is
## Pearson r across edges. Run for every FRAME lambda_m fit found, on edges
## touching sampled nodes (primary) and on all edges.
##
## R libraries: windscape 2.0.0 and terra 1.9-50 (CRAN) in BASE_DIR/software/R_libs_windscape,
##   installed with install.packages(c("windscape", "terra"), lib = <that dir>)
## Usage: Rscript frame_wind_asymmetry.R
## Output (in migration/FEEMS/): frame_wind_asymmetry_results.tsv,
##   frame_wind_asymmetry_null.tsv, frame_wind_asymmetry.jpg

BASE_DIR <- "/private/groups/meyerlab/esteiger"
.libPaths(c(file.path(BASE_DIR, "software/R_libs_windscape"),
            file.path(BASE_DIR, "software/R_libs/")))
library(windscape)
library(terra)

FEEMS_DIR <- file.path(BASE_DIR, "migration/FEEMS")
WS_DIR <- file.path(BASE_DIR, "wind/windscape")
PREFIX <- "strauss_all_maf"
N_ROT <- 200

# ---------- 1. Grid, edges, FRAME fits ----------

grid <- as.matrix(read.table(file.path(FEEMS_DIR, "grid.csv")))
colnames(grid) <- c("lon", "lat")
edges <- as.matrix(read.table(file.path(FEEMS_DIR, "edges.csv"))) - 1   # edges.csv is 1-based; FRAME tables are 0-based
edges <- cbind(a = pmin(edges[, 1], edges[, 2]), b = pmax(edges[, 1], edges[, 2]))
stopifnot(min(edges) == 0, max(edges) == nrow(grid) - 1)

fit_files <- list.files(FEEMS_DIR, pattern = sprintf("^%s_frame_lambda.*_edges\\.tsv$", PREFIX), full.names = TRUE)
lambdas <- sub(sprintf("^%s_frame_lambda(.*)_edges\\.tsv$", PREFIX), "\\1", basename(fit_files))
fit_files <- fit_files[order(as.numeric(lambdas))]; lambdas <- sort(as.numeric(lambdas))
cat("FRAME fits:", paste(lambdas, collapse = ", "), "\n")

geneflow_asym <- function(file) {
  e <- read.delim(file)
  gf <- setNames(e$m, paste(e$geneflow_from, e$geneflow_to))
  log(gf[paste(edges[, "a"], edges[, "b"])] / gf[paste(edges[, "b"], edges[, "a"])])
}
gf_asym <- sapply(fit_files, geneflow_asym)
colnames(gf_asym) <- paste0("lambda_", lambdas)
stopifnot(!anyNA(gf_asym))

nodes <- read.delim(sub("_edges\\.tsv$", "_nodes.tsv", fit_files[1]))
observed <- nodes$node[nodes$n_haps > 0]
near <- edges[, "a"] %in% observed | edges[, "b"] %in% observed
cat(sprintf("%d grid nodes with plants; %d of %d edges touch one\n", length(observed), sum(near), nrow(edges)))

# ---------- 2. Wind asymmetry: observed field and rotated fields ----------

series <- rast(file.path(WS_DIR, "windninja_tranquillon_wind_series.tif"))
n_t <- nlyr(series) / 2
u <- series[[1:n_t]]; v <- series[[(n_t + 1):(2 * n_t)]]

wind_asym <- function(theta_deg) {
  th <- theta_deg * pi / 180
  rot <- wind_series(c(u * cos(th) - v * sin(th), u * sin(th) + v * cos(th)))
  hrs <- pairwise_least_cost(wind_rose(rot, trans = 1), grid)
  log(hrs[cbind(edges[, "b"] + 1, edges[, "a"] + 1)] / hrs[cbind(edges[, "a"] + 1, edges[, "b"] + 1)])
}

# Rotations depend only on the wind and the grid, so they are cached.
# Delete the cache file to recompute.
null_file <- file.path(FEEMS_DIR, "frame_wind_asymmetry_null.tsv")
if (file.exists(null_file)) {
  cached <- read.delim(null_file, check.names = FALSE)
  stopifnot(all(cached$a == edges[, "a"]), all(cached$b == edges[, "b"]))
  w_asym <- as.matrix(cached[, grep("^rot_", names(cached))])
  cat(sprintf("Read %d cached wind rotations from %s\n", ncol(w_asym), null_file))
} else {
  angles <- seq(0, 360, length.out = N_ROT + 1)[1:N_ROT]   # angle 0 = observed wind
  t0 <- Sys.time()
  w_asym <- sapply(seq_along(angles), function(k) {
    if (k %% 20 == 0) cat(sprintf("  rotation %d / %d [%.0fs]\n", k, N_ROT, as.numeric(Sys.time() - t0, units = "secs")))
    wind_asym(angles[k])
  })
  colnames(w_asym) <- paste0("rot_", round(angles, 1))
  write.table(cbind(edges, near = near, round(w_asym, 5)), null_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---------- 3. Statistics ----------

results <- do.call(rbind, lapply(colnames(gf_asym), function(lam) {
  do.call(rbind, lapply(c("near_samples", "all_edges"), function(set) {
    keep <- if (set == "near_samples") near else rep(TRUE, nrow(edges))
    r <- apply(w_asym[keep, , drop = FALSE], 2, function(w) cor(gf_asym[keep, lam], w))
    data.frame(lambda_m = sub("lambda_", "", lam), edge_set = set, n_edges = sum(keep),
               r = r[1], p_one_sided = mean(r >= r[1]), p_two_sided = mean(abs(r) >= abs(r[1])),
               null_mean = mean(r[-1]), null_sd = sd(r[-1]))
  }))
}))
rownames(results) <- NULL
write.table(results, file.path(FEEMS_DIR, "frame_wind_asymmetry_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(results, digits = 3)

# ---------- 4. Figure: observed scatter and rotation null per lambda ----------

jpeg(file.path(FEEMS_DIR, "frame_wind_asymmetry.jpg"), width = 3.2 * length(lambdas), height = 6.2, units = "in", res = 300)
par(mfrow = c(2, length(lambdas)), mar = c(4.2, 4.2, 2.2, 0.8))
for (lam in colnames(gf_asym)) {
  plot(w_asym[near, 1], gf_asym[near, lam], pch = 19, col = adjustcolor("#009e73", 0.6), cex = 0.7,
       xlab = "Wind-flow asymmetry, ln(a→b / b→a)", ylab = "Gene-flow asymmetry, ln(a→b / b→a)",
       main = sprintf("FRAME lambda_m = %s", sub("lambda_", "", lam)), cex.main = 0.9)
  abline(h = 0, v = 0, col = "grey70", lty = 2)
  abline(lm(gf_asym[near, lam] ~ w_asym[near, 1]), lwd = 2)
}
for (lam in colnames(gf_asym)) {
  r <- apply(w_asym[near, , drop = FALSE], 2, function(w) cor(gf_asym[near, lam], w))
  hist(r[-1], breaks = 25, col = "grey85", border = "white", xlim = range(r),
       main = sprintf("r = %.3f, one-sided p = %.3f", r[1], mean(r >= r[1])), cex.main = 0.9,
       xlab = "r under rotated wind")
  abline(v = r[1], col = "#d9620b", lwd = 2)
}
dev.off()
cat("Saved results and figure to", FEEMS_DIR, "\n")
