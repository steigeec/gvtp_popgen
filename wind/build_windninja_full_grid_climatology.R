## WindNinja produces a full downscaled speed+direction GRID at every
## timestep (not just point values), but build_hybrid_windrose.py only ever
## sampled that grid at the 6 population coordinates (all that the AEM+RDA
## analysis needed). This script instead averages the WindNinja output
## across all 1464 timesteps CELL-BY-CELL over the whole domain, producing a
## genuine spatially-continuous WindNinja wind climatology -- comparable to
## hrrr_windrose_summary.txt's per-cell climatology -- so the map can show
## the full downscaled field instead of only 6 isolated corrected points.

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(terra)

WN_DIR <- "/private/groups/meyerlab/esteiger/wind/windninja_full/grib2"
OUT_DIR <- "/private/groups/meyerlab/esteiger/wind/windninja_full"

vel_files <- sort(list.files(WN_DIR, pattern = "_vel\\.asc$", recursive = TRUE, full.names = TRUE))
ang_files <- sub("_vel\\.asc$", "_ang.asc", vel_files)
stopifnot(all(file.exists(ang_files)))
cat(sprintf("Found %d WindNinja timestep grids\n", length(vel_files)))

# All grids share the same 460 m mesh (same DEM, same domain, same
# mesh_resolution every run), so they stack directly without resampling.
ref <- rast(vel_files[1])
crs(ref) <- "EPSG:32610"  # WindNinja ascii output has no embedded CRS; it's the UTM 10N mesh from the DEM

sum_u <- rast(ref); values(sum_u) <- 0
sum_v <- rast(ref); values(sum_v) <- 0
sum_speed <- rast(ref); values(sum_speed) <- 0
n <- length(vel_files)

t0 <- Sys.time()
for (i in seq_len(n)) {
  speed <- rast(vel_files[i])
  dir_from <- rast(ang_files[i])
  crs(speed) <- crs(sum_u); crs(dir_from) <- crs(sum_u)
  u <- -speed * sin(dir_from * pi / 180)
  v <- -speed * cos(dir_from * pi / 180)
  sum_u <- sum_u + u
  sum_v <- sum_v + v
  sum_speed <- sum_speed + speed
  if (i %% 200 == 0) cat(sprintf("  %d / %d [%.0fs]\n", i, n, as.numeric(Sys.time() - t0, units = "secs")))
}

mean_u <- sum_u / n
mean_v <- sum_v / n
mean_speed <- sum_speed / n
resultant_speed <- sqrt(mean_u^2 + mean_v^2)
# NOTE: terra's SpatRaster %% does not wrap negatives into [0,360) the way
# base R's %% does on plain numeric vectors (confirmed empirically) -- use
# ifel() instead of a bare %% 360 here.
resultant_dir_raw <- atan2(-mean_u, -mean_v) * 180 / pi
resultant_dir <- ifel(resultant_dir_raw < 0, resultant_dir_raw + 360, resultant_dir_raw)
R <- resultant_speed / mean_speed

out <- c(mean_speed, resultant_dir, R)
names(out) <- c("mean_speed_ms", "resultant_dir_deg_from", "resultant_R")
writeRaster(out, file.path(OUT_DIR, "windninja_full_grid_climatology.tif"), overwrite = TRUE)
cat("\nSaved ->", file.path(OUT_DIR, "windninja_full_grid_climatology.tif"), "\n")
