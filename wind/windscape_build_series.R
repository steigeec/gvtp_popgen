## Build the bloom-season u/v wind time series used by windscape_connectivity.R:
## raw HRRR, WindNinja over the region, and a finer WindNinja run over Tranquillon Mtn.
##
## Inputs:
##   hrrr/hrrr_wind_raw_points_earth.csv (download_hrrr_wind.py, rotated to
##     earth-relative by rotate_hrrr_wind.py); 478 cells x 4896 timesteps
##     (2015-2022, Jun-Oct, 00/06/12/18 UTC)
##   windninja_landfire/<domain>/*_vel.asc, *_ang.asc (run_windninja_landfire_batch.sh;
##     UTM 10N; speed in m/s, direction the wind blows from, relative to true north)
## Outputs (in BASE_DIR/wind/windscape/):
##   hrrr_wind_series.tif, windninja_wind_series.tif,
##   windninja_tranquillon_wind_series.tif
##   Each has 2 x 4896 layers ordered "uuvv" (as windscape expects), named u_/v_<UTC timestamp>.
## Grids: HRRR 0.03 deg; WindNinja region 460 m mesh to 0.005 deg; Tranquillon
##   100 m mesh to 0.001 deg, inset ~1 km from the WindNinja domain edge.
## R libraries: R_libs_windscape, R_libs (terra, data.table)
## Usage: Rscript windscape_build_series.R

BASE_DIR <- "/private/groups/meyerlab/esteiger"
.libPaths(c(file.path(BASE_DIR, "software/R_libs_windscape"),
            file.path(BASE_DIR, "software/R_libs/")))
library(terra)
library(data.table)

WIND_DIR <- file.path(BASE_DIR, "wind")
OUT_DIR <- file.path(WIND_DIR, "windscape")
dir.create(OUT_DIR, showWarnings = FALSE)

# HRRR native Lambert Conformal grid (spherical earth, 3 km spacing).
HRRR_LCC <- "+proj=lcc +lat_0=38.5 +lon_0=-97.5 +lat_1=38.5 +lat_2=38.5 +x_0=0 +y_0=0 +R=6371229 +units=m +no_defs"
HRRR_DX <- 3000

# Output grids sit inside each source domain so edge cells are not extrapolated.
HRRR_EXT <- ext(-120.72, -120.09, 34.37, 34.94)
HRRR_RES <- 0.03
WN_DOMAINS <- list(
  region      = list(ext = ext(-120.735, -120.115, 34.41, 34.94), res = 0.005, mesh = 460,
                     out = "windninja_wind_series.tif"),
  tranquillon = list(ext = ext(-120.558, -120.496, 34.559, 34.609), res = 0.001, mesh = 100,
                     out = "windninja_tranquillon_wind_series.tif")
)

# ---------- 1. HRRR points onto the Lambert grid ----------

h <- fread(file.path(WIND_DIR, "hrrr", "hrrr_wind_raw_points_earth.csv"))
h[, lon := lon - 360]
times <- sort(unique(h$datetime))
n_t <- length(times)
cat(sprintf("HRRR: %d rows, %d timesteps\n", nrow(h), n_t))

cells <- unique(h[, .(lat, lon)])
xy <- project(as.matrix(cells[, .(lon, lat)]), from = "EPSG:4326", to = HRRR_LCC)
col_f <- (xy[, 1] - min(xy[, 1])) / HRRR_DX
row_f <- (max(xy[, 2]) - xy[, 2]) / HRRR_DX
snap_err <- max(abs(c(col_f - round(col_f), row_f - round(row_f)))) * HRRR_DX
cat(sprintf("HRRR cells snap to the 3 km Lambert grid within %.1f m\n", snap_err))
stopifnot(snap_err < 100)

n_row <- max(round(row_f)) + 1
n_col <- max(round(col_f)) + 1
x0 <- min(xy[, 1]) - HRRR_DX / 2
y1 <- max(xy[, 2]) + HRRR_DX / 2
lcc <- rast(nrows = n_row, ncols = n_col, crs = HRRR_LCC,
            extent = ext(x0, x0 + n_col * HRRR_DX, y1 - n_row * HRRR_DX, y1))
cells[, cell := cellFromRowCol(lcc, round(row_f) + 1, round(col_f) + 1)]
h <- merge(h, cells, by = c("lat", "lon"))

to_layers <- function(var) {
  w <- dcast(h, cell ~ datetime, value.var = var)
  m <- matrix(NA_real_, ncell(lcc), n_t)
  m[w$cell, ] <- as.matrix(w[, as.character(times), with = FALSE])
  r <- rast(lcc, nlyrs = n_t)
  values(r) <- m
  r
}
u_lcc <- to_layers("u10")
v_lcc <- to_layers("v10")

stamp <- format(as.POSIXct(times, tz = "UTC"), "%Y%m%dT%H")

# ---------- 2. HRRR series on the lon/lat grid ----------

hrrr_tmpl <- rast(HRRR_EXT, resolution = HRRR_RES, crs = "EPSG:4326")
hrrr_series <- project(c(u_lcc, v_lcc), hrrr_tmpl, method = "bilinear")
names(hrrr_series) <- c(paste0("u_", stamp), paste0("v_", stamp))
cat(sprintf("HRRR series: %d x %d cells, %d layers, NA cells: %d\n",
            nrow(hrrr_series), ncol(hrrr_series), nlyr(hrrr_series),
            sum(is.na(values(hrrr_series[[1]])))))
writeRaster(hrrr_series, file.path(OUT_DIR, "hrrr_wind_series.tif"), overwrite = TRUE)

# ---------- 3. WindNinja series ----------

build_wn_series <- function(domain) {
  d <- WN_DOMAINS[[domain]]
  tmpl <- rast(d$ext, resolution = d$res, crs = "EPSG:4326")
  vel_files <- list.files(file.path(WIND_DIR, "windninja_landfire", domain),
                          pattern = "_vel\\.asc$", full.names = TRUE)
  # File names carry local (Pacific) time, e.g.
  # region_landfire_05-31-2022_1700_460m_non_neutral_stability_vel.asc
  m <- regmatches(basename(vel_files), regexec(
    sprintf("_(\\d{2})-(\\d{2})-(\\d{4})_(\\d{2})(\\d{2})_%dm_non_neutral_stability_vel\\.asc$", d$mesh),
    basename(vel_files)))
  local <- as.POSIXct(sapply(m, function(x) sprintf("%s-%s-%s %s:%s", x[4], x[2], x[3], x[5], x[6])),
                      format = "%Y-%m-%d %H:%M", tz = "America/Los_Angeles")
  vel_by_stamp <- setNames(vel_files, format(local, "%Y%m%dT%H", tz = "UTC"))
  stopifnot(all(stamp %in% names(vel_by_stamp)))
  cat(sprintf("WindNinja %s: %d grids, all %d HRRR timesteps matched\n", domain, length(vel_files), n_t))

  u_wn <- v_wn <- vector("list", n_t)
  t0 <- Sys.time()
  for (k in seq_len(n_t)) {
    vel <- rast(vel_by_stamp[stamp[k]])
    ang <- rast(sub("_vel\\.asc$", "_ang.asc", vel_by_stamp[stamp[k]]))
    rad <- ang * pi / 180
    uv <- c(-vel * sin(rad), -vel * cos(rad))  # ang is direction from, so negate
    if (is.na(crs(uv)) || crs(uv) == "") crs(uv) <- "EPSG:32610"
    uv <- project(uv, tmpl, method = "bilinear")
    u_wn[[k]] <- uv[[1]]
    v_wn[[k]] <- uv[[2]]
    if (k %% 300 == 0) cat(sprintf("  %d / %d [%.0fs]\n", k, n_t, as.numeric(Sys.time() - t0, units = "secs")))
  }
  series <- c(rast(u_wn), rast(v_wn))
  names(series) <- c(paste0("u_", stamp), paste0("v_", stamp))
  cat(sprintf("WindNinja %s series: %d x %d cells, %d layers, NA cells: %d\n", domain,
              nrow(series), ncol(series), nlyr(series), sum(is.na(values(series[[1]])))))
  writeRaster(series, file.path(OUT_DIR, d$out), overwrite = TRUE)
}

for (domain in names(WN_DOMAINS)) build_wn_series(domain)

cat("Saved ->", OUT_DIR, "\n")
