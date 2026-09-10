## Wind-mediated dispersal / gene-flow test -- INTERIM, non-directional version.
## Companion to structure/mantel_IBD.R (plain isolation-by-distance).
##
## Hypothesis: wind-mediated pollen/seed dispersal shapes population
## differentiation independently of (or beyond) simple geographic distance,
## because Point Conception's coastal bend produces a sharp wind
## acceleration/direction gradient (confirmed empirically below) rather than
## a smooth distance-decay field.
##
## IMPORTANT LIMITATION: Fst is symmetric by construction (built from
## allele-frequency variance; it cannot distinguish "gene flow went A->B"
## from "gene flow went B->A"). wind_connectivity_sym below is a symmetrized
## downwind transport score, tested against symmetric Fst -- so this script
## can only test TOTAL wind connectivity between a pair, not directionality.
## A genuinely directional test needs either (a) an asymmetric genetic
## differentiation estimate (e.g. divMigrate-style relative migration --
## caveat: Gst-asymmetry can reflect asymmetric drift/Ne, not just asymmetric
## migration), or (b) AEM + RDA on population allele frequencies with
## directional wind-weighted spatial eigenvectors (Blanchet et al. 2008) --
## the textbook method for directional processes, no drift/Ne confound, but
## a different response variable (allele frequencies, not pairwise Fst).
## This has since been built separately: see wind_aem_rda.R.
##
## BLOOM-SEASON RESTRICTION: wind only plausibly matters here via pollen/seed
## transport, which can only happen while the plant is flowering (June-
## September). The station wind-rose statistics feeding wind_support/
## wind_connectivity_sym (windrose_summary.txt, windrose_sectors.txt; see
## build_windrose_table.py) are therefore computed on Jun-Sep observations
## only, not an annual average -- this matters empirically, not just in
## principle: restricting to bloom season raised the directional steadiness
## (resultant_R) at every station (e.g. PTGC1 0.80 -> 0.93, TQS 0.49 -> 0.75),
## consistent with California's coastal NW wind regime being strongest and
## most persistent in summer.
##
## RESOLVED GAP: the GWA wind-speed layer (download_gwa_windspeed.py) is a
## long-term *annual* average with no seasonal breakdown available -- HRRR
## replaces it for both the map and this analysis with a bloom-season-only,
## spatially continuous speed+direction field (3 km, 478 grid cells). GWA's
## 250 m resolution is still finer than HRRR's 3 km, so GWA remains the
## better source if a purely-speed, annual-average exposure surface is
## wanted for some other purpose -- it's just no longer used here.
##
## wind_connectivity_sym[i,j] = symmetrized downwind transport potential
##     between i and j, built from each site's own local wind rose (direction
##     x speed frequency table) and the compass bearing between the two
##     sites: wind_connectivity_sym = (wind_support(i,j) + wind_support(j,i)) / 2.
##     The raw asymmetric matrix is still saved (wind_support_asymmetric.txt)
##     as the input for whichever directional method is built next.
##
## UPDATE: wind rose now comes from HRRR (3 km, hourly, 2022-2024 Jun-Sep --
## see download_hrrr_wind.py / process_hrrr_climatology.py), not the 4
## weather stations. This directly fixes the previous limitation of only
## having real direction data at 4 point locations (nearest station could be
## 1.5-32.7 km from a population, e.g. CDFW's nearest was SBA at 32.7 km): with
## HRRR, each of the 6 populations gets its OWN nearest grid cell out of 478
## covering the whole study area, so no population borrows direction data
## from a distant point anymore. See station_map history below for what this
## replaced; the old station-based windrose_summary.txt/windrose_sectors.txt
## (build_windrose_table.py) are left as-is for reference/comparison but are
## no longer read by this script.
##
## [historical] Site -> nearest wind STATION assignment (superseded above):
##   VAFB_Arguello, Dangermond, Point_Conception -> PTGC1 (Pt. Arguello C-MAN,
##     1.5-18.7 km away; PTGC1 is literally at the point)
##   VAFB_Lion_Head -> VGN (Vandenberg SFB North ASOS, 1.6 km)
##   Strauss -> TQS (Vandenberg SFB South ASOS, 10.6 km)
##   CDFW -> SBA (Santa Barbara Muni ASOS, 32.7 km -- best available; CDFW
##     sits in the lee of the cape and no closer station samples that
##     sheltered regime. Known limitation, not a close match.)

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(vegan)
library(geosphere)
library(terra)

WIND_DIR <- "/private/groups/meyerlab/esteiger/wind"

# ---------- 1. Site centroids (same construction as structure/mantel_IBD.R) ----------

field <- read.delim("/private/groups/meyerlab/esteiger/field_data.txt", header = TRUE, sep = "\t")
centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field, FUN = mean)

name_map <- c(
  "Point_Conception" = "conception",
  "CDFW"             = "cdfw",
  "VAFB_Arguello"    = "arguello",
  "Strauss"          = "strauss",
  "VAFB_Lion_Head"   = "lionhead",
  "Dangermond"       = "Dangermond"
)
centroids$pop <- name_map[centroids$Clusters]

# ---------- 2. Pairwise Fst (nuclear, population-level, weighted) ----------

fst_df <- read.delim(
  "/private/groups/meyerlab/esteiger/structure/Fst/populations/pairwise_auto/pairwise_Fst_populations.txt",
  header = TRUE, sep = "\t"
)
pops <- sort(unique(c(fst_df$pop1, fst_df$pop2)))
n <- length(pops)
stopifnot(nrow(fst_df) == n * (n - 1) / 2)

fst_mat <- matrix(0, n, n, dimnames = list(pops, pops))
for (i in seq_len(nrow(fst_df))) {
  fst_mat[fst_df$pop1[i], fst_df$pop2[i]] <- fst_df$Fst_weighted[i]
  fst_mat[fst_df$pop2[i], fst_df$pop1[i]] <- fst_df$Fst_weighted[i]
}
fst_lin <- fst_mat / (1 - fst_mat)  # Rousset 1997 linearization

centroids <- centroids[match(pops, centroids$pop), ]
rownames(centroids) <- centroids$pop

# ---------- 3. Geographic distance (km) + log distance ----------

coords <- as.matrix(centroids[, c("Long", "Lat")])
rownames(coords) <- centroids$pop
geo_mat <- as.matrix(distm(coords, fun = distHaversine)) / 1000
dimnames(geo_mat) <- list(pops, pops)
geo_log <- log(geo_mat)
diag(geo_log) <- 0

# ---------- 4. Directional wind connectivity from HRRR grid-cell wind roses ----------

hrrr_summary <- read.delim(file.path(WIND_DIR, "hrrr", "hrrr_windrose_summary.txt"), header = TRUE, sep = "\t")
sectors <- read.delim(file.path(WIND_DIR, "hrrr", "hrrr_windrose_sectors.txt"), header = TRUE, sep = "\t")

# Each population gets its OWN nearest HRRR grid cell (out of 478 covering
# the whole study area) rather than borrowing direction data from one of
# only 4 distant stations.
hrrr_summary$lon180 <- ((hrrr_summary$lon + 180) %% 360) - 180  # HRRR 0-360 -> -180/180
hrrr_lonlat <- as.matrix(hrrr_summary[, c("lon180", "lat")])
nearest_cell <- function(lon, lat) {
  d <- distHaversine(matrix(c(lon, lat), ncol = 2), hrrr_lonlat)
  list(cell_id = hrrr_summary$cell_id[which.min(d)], dist_km = min(d) / 1000)
}
nc <- mapply(nearest_cell, centroids$Long, centroids$Lat, SIMPLIFY = FALSE)
centroids$cell_id <- vapply(nc, `[[`, character(1), "cell_id")
cat("Population -> nearest HRRR grid cell (distance, km):\n")
for (p in seq_len(nrow(centroids))) {
  cat(sprintf("  %-20s -> %s (%.2f km)\n", centroids$pop[p], nc[[p]]$cell_id, nc[[p]]$dist_km))
}

# Half-cosine directional kernel: a sector contributes to "wind support"
# toward bearing b in proportion to cos(angle between them) when that angle
# is < 90 deg, and not at all beyond +/-90 deg (wind blowing away from the
# target contributes nothing). Weights are normalized to sum to 1 so the
# result is a frequency-and-speed-weighted score in the same units as
# mean_speed_ms, interpretable as "effective downwind transport speed
# toward that bearing".
wind_support <- function(cell_id, bearing_deg) {
  sec <- sectors[sectors$cell_id == cell_id, ]
  ang <- abs((sec$sector_center_deg - bearing_deg + 180) %% 360 - 180)
  w <- pmax(0, cos(ang * pi / 180))
  if (sum(w) == 0) return(0)
  w <- w / sum(w)
  sum(w * sec$frequency * sec$mean_speed_ms)
}

wind_support_asym <- matrix(0, n, n, dimnames = list(pops, pops))
for (i in pops) {
  for (j in pops) {
    if (i == j) next
    b <- bearing(coords[i, ], coords[j, ]) %% 360
    wind_support_asym[i, j] <- wind_support(centroids[i, "cell_id"], b)
  }
}
write.table(round(wind_support_asym, 4), file.path(WIND_DIR, "wind_support_asymmetric.txt"),
            sep = "\t", quote = FALSE)

wind_connectivity_sym <- (wind_support_asym + t(wind_support_asym)) / 2

cat("Asymmetry check -- |wind_support(i,j) - wind_support(j,i)|, upper triangle:\n")
asym_diff <- abs(wind_support_asym - t(wind_support_asym))
print(round(asym_diff[upper.tri(asym_diff)], 3))
cat("(This asymmetry is exactly what this script CANNOT test against symmetric Fst --\n")
cat(" see header. Reported here only as a diagnostic of how much directional\n")
cat(" information the symmetrization below is discarding.)\n")

# ---------- 5. Mantel tests ----------

fst_dist <- as.dist(fst_lin)
geo_dist <- as.dist(geo_log)
conn_dist <- as.dist(-wind_connectivity_sym)  # negate: higher connectivity -> treated as "closer"

set.seed(42)
cat("\n=== Mantel: Fst ~ geographic distance (reference; matches structure/mantel_IBD.R) ===\n")
m_geo <- mantel(fst_dist, geo_dist, method = "pearson", permutations = 9999)
cat(sprintf("r = %.4f, p = %.4f\n", m_geo$statistic, m_geo$signif))

cat("\n=== Partial Mantel: Fst ~ wind_connectivity_sym | geographic distance ===\n")
pm_conn <- mantel.partial(fst_dist, conn_dist, geo_dist, method = "pearson", permutations = 9999)
cat(sprintf("r = %.4f, p = %.4f  (conn_dist = -connectivity, so a POSITIVE r here means\n",
            pm_conn$statistic, pm_conn$signif))
cat("  higher wind connectivity -> lower Fst, as predicted)\n")

# ---------- 6. MMRR: Fst ~ geo_dist + wind_connectivity_sym ----------
# Self-contained implementation (Wang 2013 permutation scheme: permute the
# row/column order of the response matrix, refit, build a null distribution
# of coefficients and R^2). Avoids adding a new package dependency.

mmrr <- function(Y, Xlist, nperm = 9999, seed = 42) {
  n <- nrow(Y)
  lower <- lower.tri(Y)
  y <- Y[lower]
  X <- sapply(Xlist, function(m) m[lower])
  X <- cbind(Intercept = 1, X)
  fit <- lm.fit(X, y)
  coef_obs <- fit$coefficients
  r2 <- function(yv, yhat) 1 - sum((yv - yhat)^2) / sum((yv - mean(yv))^2)
  r2_obs <- r2(y, X %*% coef_obs)

  set.seed(seed)
  perm_coef <- matrix(NA, nperm, length(coef_obs))
  perm_r2 <- numeric(nperm)
  for (p in seq_len(nperm)) {
    idx <- sample(n)
    yp <- Y[idx, idx][lower]
    fit_p <- lm.fit(X, yp)
    perm_coef[p, ] <- fit_p$coefficients
    perm_r2[p] <- r2(yp, X %*% fit_p$coefficients)
  }
  pvals <- vapply(seq_along(coef_obs), function(i) mean(abs(perm_coef[, i]) >= abs(coef_obs[i])), numeric(1))
  list(coefficients = coef_obs, pvalues = setNames(pvals, names(coef_obs)),
       r2 = r2_obs, r2_pvalue = mean(perm_r2 >= r2_obs), nperm = nperm)
}

cat("\n=== MMRR: Fst ~ geo_dist + wind_connectivity_sym ===\n")
mm <- mmrr(fst_lin, list(geo_dist = geo_log, wind_connectivity_sym = -wind_connectivity_sym), nperm = 9999)
cat(sprintf("R^2 = %.4f (p = %.4f)\n", mm$r2, mm$r2_pvalue))
for (nm in names(mm$coefficients)) {
  cat(sprintf("  %-25s coef = %8.5f  p = %.4f\n", nm, mm$coefficients[nm], mm$pvalues[nm]))
}

# ---------- 7. Save summary + diagnostic plot ----------

summary_txt <- file.path(WIND_DIR, "wind_connectivity_results.txt")
sink(summary_txt)
cat("Wind-mediated gene flow test -- summary (interim, non-directional; see script header)\n")
cat("========================================================================\n\n")
cat(sprintf("Mantel Fst ~ geo dist:                          r = %.4f, p = %.4f\n", m_geo$statistic, m_geo$signif))
cat(sprintf("Partial Mantel Fst ~ wind connectivity | geo:   r = %.4f, p = %.4f\n", pm_conn$statistic, pm_conn$signif))
cat(sprintf("\nMMRR R^2 = %.4f (p = %.4f)\n", mm$r2, mm$r2_pvalue))
for (nm in names(mm$coefficients)) {
  cat(sprintf("  %-25s coef = %8.5f  p = %.4f\n", nm, mm$coefficients[nm], mm$pvalues[nm]))
}
sink()
cat("\nSaved summary ->", summary_txt, "\n")

jpeg(file.path(WIND_DIR, "wind_connectivity_plot.jpg"), width = 5, height = 5, units = "in", res = 300)
plot(as.vector(-conn_dist), as.vector(fst_dist),
     xlab = "Symmetrized wind connectivity", ylab = expression(F[ST]/(1-F[ST])),
     pch = 21, bg = "#d55e00", main = sprintf("partial r=%.2f, p=%.3f", pm_conn$statistic, pm_conn$signif))
dev.off()
cat("Saved plot -> wind_connectivity_plot.jpg\n")

# ---------- 9. Wind map ----------
# HRRR bloom-season (Jun-Sep, 2022-2024) speed + direction, spatially
# continuous across all 478 grid cells covering the study area -- replaces
# the earlier GWA-speed-raster + 4-station-arrow map (GWA is an annual
# average with no direction; stations gave direction at only 4 points, up to
# 32.7 km from a population). HRRR is coarser (3 km vs. GWA's 250 m) but is
# bloom-season-specific and self-consistent (speed and direction from the
# same model grid, not two different sources).

bbox_extent <- ext(-120.7, -120.1, 34.4, 34.9)

hrrr_pts <- hrrr_summary
hrrr_pts$lon <- hrrr_pts$lon180  # already converted to -180/180 in section 4
pts_v <- vect(hrrr_pts, geom = c("lon", "lat"), crs = "EPSG:4326")

# Template grid at ~HRRR's native spacing (3 km) so each grid cell maps to
# ~one raster cell (rather than a much finer template that would leave most
# cells empty, or a much coarser one that would average multiple HRRR cells
# together).
template <- rast(bbox_extent, resolution = 0.03, crs = "EPSG:4326")
speed_r <- rasterize(pts_v, template, field = "mean_speed_ms", fun = "mean")
# HRRR's native grid doesn't align exactly with this lat/lon template, so a
# few template cells get no point and come out NA ("holes" in the map). Fill
# those with the mean of their non-NA neighbors (one pass is enough given
# how sparse the gaps are relative to cell size).
speed_r <- focal(speed_r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)

# Arrow length in degrees per m/s -- chosen so a ~3-8 m/s wind draws a
# visually reasonable arrow relative to this bbox's ~0.5-0.6 deg span.
SPEED_SCALE <- 0.006
# "resultant_dir_deg_from" is the compass bearing the wind blows FROM (met.
# convention); arrows are drawn pointing where the wind blows TO, i.e.
# bearing + 180.
to_bearing_rad <- (hrrr_pts$resultant_dir_deg_from + 180) %% 360 * pi / 180
hrrr_pts$dx <- sin(to_bearing_rad) * hrrr_pts$mean_speed_ms * SPEED_SCALE
hrrr_pts$dy <- cos(to_bearing_rad) * hrrr_pts$mean_speed_ms * SPEED_SCALE

# Subsample the 478 cells for arrow drawing (every other grid cell in each
# direction) so arrows are legible rather than overlapping at 3 km spacing.
in_bbox <- hrrr_pts$lon >= bbox_extent[1] & hrrr_pts$lon <= bbox_extent[2] &
           hrrr_pts$lat >= bbox_extent[3] & hrrr_pts$lat <= bbox_extent[4]
hp <- hrrr_pts[in_bbox, ]
keep <- (round((hp$lon - min(hp$lon)) / 0.03) %% 2 == 0) &
        (round((hp$lat - min(hp$lat)) / 0.03) %% 2 == 0)
hp <- hp[keep, ]

# Display labels for the map (place names), keyed by the internal pop codes
# used throughout this script.
pop_labels <- c(
  lionhead   = "Lion's Head",
  arguello   = "Pt. Arguello",
  strauss    = "Tranquillon",
  Dangermond = "Dangermond",
  conception = "Pt. Conception",
  cdfw       = "Gaviota"
)

jpeg(file.path(WIND_DIR, "wind_map.jpg"), width = 7, height = 6.3, units = "in", res = 300)
# Palette truncated to its lighter ~75%: dropping the darkest end keeps the
# black wind arrows visible against the raster instead of disappearing into
# near-black cells at high wind speed.
speed_col <- hcl.colors(140, "YlGnBu", rev = TRUE)[1:100]
# xlab/ylab left blank and added via mtext() below -- terra's plot() manages
# its own layout (it resets par(mar) to reserve space for the color-scale
# legend, overriding anything set beforehand), so axis titles are added
# afterward with line= values measured to fit the margin terra actually left.
plot(speed_r, ext = bbox_extent,
     col = speed_col,
     xlab = "", ylab = "",
     cex.axis = 1.3)
mtext("Longitude", side = 1, line = 2.1, cex = 1.3)
mtext("Latitude", side = 2, line = 2.2, cex = 1.3)

# Wind-rose arrows: length ~ mean speed, opacity ~ resultant_R (directional
# steadiness) -- steadier prevailing wind draws a more solid arrow.
arrow_col <- rgb(0, 0, 0, alpha = pmax(0.3, hp$resultant_R))
arrows(hp$lon, hp$lat, hp$lon + hp$dx, hp$lat + hp$dy,
       length = 0.05, lwd = 1.2, col = arrow_col)

# Population sites. Arguello and Tranquillon (Strauss) sit close enough
# together that labeling both directly above their points collides.
# Arguello also sits near the left edge of the bbox, so its label goes
# below (pos=1) rather than to the side, which would run off the plot.
label_pos <- setNames(rep(3, nrow(coords)), rownames(coords))
label_pos["arguello"] <- 1

points(coords[, "Long"], coords[, "Lat"], pch = 19, col = "black", cex = 1.4)
text(coords[, "Long"], coords[, "Lat"], labels = pop_labels[rownames(coords)],
     pos = label_pos[rownames(coords)], cex = 1.1, font = 2)

dev.off()
cat("Saved wind map -> wind_map.jpg\n")
