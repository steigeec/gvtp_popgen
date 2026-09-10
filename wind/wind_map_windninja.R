## Same wind map as wind_connectivity.R section 9 (HRRR-only background field
## + coarse 3 km wind-rose arrows), with a second arrow set drawn ACROSS THE
## WHOLE DOMAIN (not just at the 6 population points -- WindNinja produces a
## full downscaled grid at every timestep, so a full-domain overlay is what
## the data actually supports) using the same hybrid convention as the
## AEM+RDA analysis: WindNinja direction (refined at 460 m, averaged
## cell-wise across all 1464 timesteps -- see
## build_windninja_full_grid_climatology.R) with HRRR speed. Drawn in bold
## black on top of the lighter original HRRR arrows so the two direction
## estimates can be compared visually everywhere, not just at 6 sites.

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(terra)

WIND_DIR <- "/private/groups/meyerlab/esteiger/wind"

# ---------- Site centroids (same construction as wind_connectivity.R) ----------

field <- read.delim("/private/groups/meyerlab/esteiger/field_data.txt", header = TRUE, sep = "\t")
centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field, FUN = mean)

name_map <- c(
  "Point_Conception" = "conception", "CDFW" = "cdfw",
  "VAFB_Arguello"    = "arguello",   "Strauss" = "strauss",
  "VAFB_Lion_Head"   = "lionhead",   "Dangermond" = "Dangermond"
)
centroids$pop <- name_map[centroids$Clusters]
rownames(centroids) <- centroids$pop
coords <- as.matrix(centroids[, c("Long", "Lat")])
rownames(coords) <- centroids$pop

# ---------- HRRR background field + coarse grid-cell arrows ----------

hrrr_summary <- read.delim(file.path(WIND_DIR, "hrrr", "hrrr_windrose_summary.txt"), header = TRUE, sep = "\t")
hrrr_summary$lon180 <- ((hrrr_summary$lon + 180) %% 360) - 180
hrrr_pts <- hrrr_summary
hrrr_pts$lon <- hrrr_pts$lon180

bbox_extent <- ext(-120.7, -120.1, 34.4, 34.9)
pts_v <- vect(hrrr_pts, geom = c("lon", "lat"), crs = "EPSG:4326")
template <- rast(bbox_extent, resolution = 0.03, crs = "EPSG:4326")
speed_r <- rasterize(pts_v, template, field = "mean_speed_ms", fun = "mean")
speed_r <- focal(speed_r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)

SPEED_SCALE <- 0.006
to_bearing_rad <- (hrrr_pts$resultant_dir_deg_from + 180) %% 360 * pi / 180
hrrr_pts$dx <- sin(to_bearing_rad) * hrrr_pts$mean_speed_ms * SPEED_SCALE
hrrr_pts$dy <- cos(to_bearing_rad) * hrrr_pts$mean_speed_ms * SPEED_SCALE

in_bbox <- hrrr_pts$lon >= bbox_extent[1] & hrrr_pts$lon <= bbox_extent[2] &
           hrrr_pts$lat >= bbox_extent[3] & hrrr_pts$lat <= bbox_extent[4]
hp <- hrrr_pts[in_bbox, ]
keep <- (round((hp$lon - min(hp$lon)) / 0.03) %% 2 == 0) &
        (round((hp$lat - min(hp$lat)) / 0.03) %% 2 == 0)
hp <- hp[keep, ]

# ---------- Hybrid (WindNinja direction + HRRR speed) arrows, FULL DOMAIN ----------
# Sample WindNinja's own full-grid direction climatology (UTM 10N) at the
# same subsampled points used for the background HRRR arrows (hp, lon/lat),
# then pair each with that point's own HRRR speed -- same hybrid convention
# used throughout (WindNinja direction + HRRR speed), just extended from 6
# population points to the whole domain.

wn_grid <- rast(file.path(WIND_DIR, "windninja_full", "windninja_full_grid_climatology.tif"))
hp_v <- vect(hp, geom = c("lon", "lat"), crs = "EPSG:4326")
hp_v_utm <- project(hp_v, crs(wn_grid))
wn_dir_at_hp <- extract(wn_grid[["resultant_dir_deg_from"]], hp_v_utm)[, 2]

hybrid_full <- hp[!is.na(wn_dir_at_hp), ]
wn_dir_valid <- wn_dir_at_hp[!is.na(wn_dir_at_hp)]

to_bearing_rad_h <- (wn_dir_valid + 180) %% 360 * pi / 180
hybrid_full$dx <- sin(to_bearing_rad_h) * hybrid_full$mean_speed_ms * SPEED_SCALE
hybrid_full$dy <- cos(to_bearing_rad_h) * hybrid_full$mean_speed_ms * SPEED_SCALE
cat(sprintf("Full-domain hybrid arrows: %d / %d subsampled points inside the WindNinja grid\n",
            nrow(hybrid_full), nrow(hp)))

# ---------- Plot ----------

pop_labels <- c(
  lionhead   = "Lion's Head", arguello = "Pt. Arguello", strauss = "Tranquillon",
  Dangermond = "Dangermond", conception = "Pt. Conception", cdfw = "Gaviota"
)
label_pos <- setNames(rep(3, nrow(coords)), rownames(coords))
label_pos["arguello"] <- 1

jpeg(file.path(WIND_DIR, "wind_map_WindNinja.jpg"), width = 7.4, height = 6.6, units = "in", res = 300)
speed_col <- hcl.colors(140, "YlGnBu", rev = TRUE)[1:100]
# Explicit mar (bottom, left, top, right): large enough for the bigger axis
# text and titles below -- terra's own auto-margin sizing (see note in
# wind_connectivity.R) is tuned for the smaller default text and clips these.
plot(speed_r, ext = bbox_extent, col = speed_col, xlab = "", ylab = "", axes = FALSE,
     mar = c(7, 8, 2.5, 7.5),
     plg = list(title = "Wind speed (m/s)", title.cex = 1.1, cex = 1.1))
# Custom axes: every-other default tick (0.2 deg spacing instead of 0.1) so
# the larger label text below has room without crowding/overlapping.
axis(1, at = seq(-120.7, -120.1, by = 0.2), cex.axis = 1.9, padj = 0.5)
axis(2, at = seq(34.4, 34.9, by = 0.2), cex.axis = 1.9, las = 1)
box()
mtext("Longitude", side = 1, line = 5, cex = 1.5)
mtext("Latitude", side = 2, line = 6, cex = 1.5)

# Original HRRR (3 km) arrows, unchanged, for context.
arrow_col <- rgb(0, 0, 0, alpha = pmax(0.3, hp$resultant_R))
arrows(hp$lon, hp$lat, hp$lon + hp$dx, hp$lat + hp$dy,
       length = 0.05, lwd = 1.2, col = arrow_col)

# WindNinja-refined (460 m) arrows across the whole domain, on top -- bold
# black (thicker, fully opaque) so they stand out from the lighter HRRR
# arrows without needing a separate legend.
arrows(hybrid_full$lon, hybrid_full$lat,
       hybrid_full$lon + hybrid_full$dx, hybrid_full$lat + hybrid_full$dy,
       length = 0.05, lwd = 1.6, col = "black")

points(coords[, "Long"], coords[, "Lat"], pch = 19, col = "black", cex = 1.4)
text(coords[, "Long"], coords[, "Lat"], labels = pop_labels[rownames(coords)],
     pos = label_pos[rownames(coords)], cex = 1.1, font = 2)

dev.off()
cat("Saved ->", file.path(WIND_DIR, "wind_map_WindNinja.jpg"), "\n")
