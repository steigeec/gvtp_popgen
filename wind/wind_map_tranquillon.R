## WindNinja-based wind map zoomed to Tranquillon Mtn (Strauss population)
## individual sample locations -- the fine-scale analog of wind_map_WindNinja.jpg
## (which covers the whole 6-population range). Unlike that map, this one
## uses WindNinja's OWN speed and direction throughout (not the HRRR-hybrid
## convention used for the AEM+RDA analysis): HRRR's 3 km grid has no real
## information at this ~2.5 km site (all individuals fall in essentially one
## HRRR cell -- see the Tranquillon microvariation check), so showing HRRR
## speed here would just be a flat, uninformative background. WindNinja's
## own 460 m field is the whole point of this map -- it resolves 18 distinct
## cells across the site (speed 3.1-5.1 m/s, direction spread ~32 deg).

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(terra)

WIND_DIR <- "/private/groups/meyerlab/esteiger/wind"

# ---------- Individual sample coordinates (Strauss / Tranquillon Mtn) ----------

field <- read.delim("/private/groups/meyerlab/esteiger/field_data.txt", header = TRUE, sep = "\t")
strauss <- field[field$Clusters == "Strauss", ]
cat(sprintf("%d individuals, lat %.4f-%.4f, lon %.4f-%.4f\n",
            nrow(strauss), min(strauss$Lat), max(strauss$Lat), min(strauss$Long), max(strauss$Long)))

# ---------- WindNinja full-grid climatology, cropped to Tranquillon + buffer ----------

wn <- rast(file.path(WIND_DIR, "windninja_full", "windninja_full_grid_climatology.tif"))

# Generous fixed buffer (~1 km, about 2 WindNinja cells) well beyond what's
# actually shown -- terra's plot() applies its own auto-padding around
# whatever ext it's given regardless of par(xaxs/yaxs) (same override
# behavior noted in wind_connectivity.R), so the crop needs to extend past
# the intended display window or that padding shows up as blank white
# instead of real data.
pad_deg <- 0.01
bbox_extent <- ext(min(strauss$Long) - pad_deg, max(strauss$Long) + pad_deg,
                    min(strauss$Lat) - pad_deg, max(strauss$Lat) + pad_deg)

bbox_v <- as.polygons(bbox_extent)
crs(bbox_v) <- "EPSG:4326"
bbox_utm <- project(bbox_v, crs(wn))
wn_crop <- crop(wn, bbox_utm)

# Reproject the cropped WindNinja grid back to lon/lat for plotting alongside
# the individual sample points (which come in lon/lat from field_data.txt).
wn_ll <- project(wn_crop, "EPSG:4326")

# Tighter display window (smaller padding than the crop above) -- real data
# now extends beyond it on every side, so terra's own auto-padding around
# this window lands on real raster cells instead of blank space.
disp_pad_lat <- 0.15 * (max(strauss$Lat) - min(strauss$Lat))
disp_pad_lon <- 0.15 * (max(strauss$Long) - min(strauss$Long))
display_extent <- ext(min(strauss$Long) - disp_pad_lon, max(strauss$Long) + disp_pad_lon,
                       min(strauss$Lat) - disp_pad_lat, max(strauss$Lat) + disp_pad_lat)

# ---------- Arrows at every WindNinja cell (native 460 m resolution -- the
#            crop is small enough that no subsampling is needed) ----------

speed_r <- wn_ll[["mean_speed_ms"]]
dir_r <- wn_ll[["resultant_dir_deg_from"]]

cell_pts <- as.points(speed_r, na.rm = TRUE)
cell_df <- data.frame(crds(cell_pts), speed = values(cell_pts)[, 1])
cell_df$dir <- extract(dir_r, cell_pts)[, 2]

SPEED_SCALE <- 0.0004  # arrow length (deg) per m/s -- tuned for this ~460m-cell-spaced grid
to_bearing_rad <- (cell_df$dir + 180) %% 360 * pi / 180
cell_df$dx <- sin(to_bearing_rad) * cell_df$speed * SPEED_SCALE
cell_df$dy <- cos(to_bearing_rad) * cell_df$speed * SPEED_SCALE

# ---------- Plot ----------

jpeg(file.path(WIND_DIR, "wind_map_tranquillon.jpg"), width = 7.4, height = 6.6, units = "in", res = 300)
speed_col <- hcl.colors(140, "YlGnBu", rev = TRUE)[1:100]
plot(speed_r, ext = display_extent, col = speed_col, xlab = "", ylab = "", axes = FALSE,
     mar = c(7, 8, 2.5, 7.5),
     plg = list(title = "Wind speed (m/s)", title.cex = 1.1, cex = 1.1))
axis(1, cex.axis = 1.7, padj = 0.5)
axis(2, cex.axis = 1.7, las = 1)
box()
mtext("Longitude", side = 1, line = 5, cex = 1.5)
mtext("Latitude", side = 2, line = 6, cex = 1.5)

arrows(cell_df[, 1], cell_df[, 2], cell_df[, 1] + cell_df$dx, cell_df[, 2] + cell_df$dy,
       length = 0.06, lwd = 2, col = "black")

# Individual sample locations (737 of the 738 field records here have
# genotype data used downstream; all 738 field coordinates are shown).
points(strauss$Long, strauss$Lat, pch = 21, bg = "white", col = "black", cex = 0.7, lwd = 0.6)

dev.off()
cat("Saved ->", file.path(WIND_DIR, "wind_map_tranquillon.jpg"), "\n")
