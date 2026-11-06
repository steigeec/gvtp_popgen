## Pairwise directional wind travel times among sites with windscape
## (Kling & Ackerly 2020, 2021), for population centroids and Tranquillon Mtn individuals.
##
## Inputs: hrrr_wind_series.tif, windninja_wind_series.tif,
##   windninja_tranquillon_wind_series.tif (windscape_build_series.R);
##   field_data.txt; strauss_individual/strauss_individuals.txt;
##   USGS 3DEP 1/3 arc-second DEM (read remotely, for land vs water)
## Outputs (in BASE_DIR/wind/windscape/):
##   <clim>_wind_rose.tif (unweighted roses, used for maps), <clim>_land_fraction.tif
##   pop_hours_hrrr.txt, pop_hours_windninja.txt (6 population centroids)
##   strauss_hours_windninja.txt (736 Tranquillon Mtn individuals, 100 m WindNinja run)
##   Matrices are hours from row site (source) to column site; asymmetric.
## R libraries: R_libs_windscape, R_libs (windscape, terra)
## Usage: Rscript windscape_connectivity.R

BASE_DIR <- "/private/groups/meyerlab/esteiger"
.libPaths(c(file.path(BASE_DIR, "software/R_libs_windscape"),
            file.path(BASE_DIR, "software/R_libs/")))
library(windscape)
library(terra)

WIND_DIR <- file.path(BASE_DIR, "wind")
WS_DIR <- file.path(WIND_DIR, "windscape")
TRANS <- 1          # conductance proportional to speed, so costs are hours
WATER_WEIGHT <- 0.1 # conductance over water down-weighted by 90%
DEM_3DEP <- "/vsicurl/https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/current/n35w121/USGS_13_n35w121.tif"

# ---------- 1. Wind roses ----------

roses <- list()
series <- list()
for (clim in c("hrrr", "windninja", "windninja_tranquillon")) {
  t0 <- Sys.time()
  series[[clim]] <- wind_series(file.path(WS_DIR, sprintf("%s_wind_series.tif", clim)))
  roses[[clim]] <- wind_rose(series[[clim]], trans = TRANS)
  writeRaster(roses[[clim]], file.path(WS_DIR, sprintf("%s_wind_rose.tif", clim)), overwrite = TRUE)
  cat(sprintf("%s wind rose: %d x %d cells [%.0fs]\n", clim, nrow(roses[[clim]]), ncol(roses[[clim]]),
              as.numeric(Sys.time() - t0, units = "secs")))
}

# ---------- 1b. Down-weight conductance over water ----------

# 3DEP codes open water as 0 m. A cell is water if less than half of it is land.
# The roses written above stay unweighted; all travel times use weighted roses.
dem <- rast(DEM_3DEP)
for (clim in names(roses)) {
  grid <- roses[[clim]][[1]]
  land <- crop(dem, ext(grid) + 0.01) > 0
  frac <- resample(land, grid, method = "average")
  names(frac) <- "land_fraction"
  writeRaster(frac, file.path(WS_DIR, sprintf("%s_land_fraction.tif", clim)), overwrite = TRUE)
  w <- ifel(frac >= 0.5, 1, WATER_WEIGHT)
  roses[[clim]] <- weight_conductance(roses[[clim]], w)
  cat(sprintf("%s: %d of %d cells water (conductance x %.1f)\n", clim,
              sum(values(w) < 1, na.rm = TRUE), ncell(w), WATER_WEIGHT))
}

# ---------- 2. Sites ----------

name_map <- c(
  "Point_Conception" = "conception", "CDFW" = "cdfw",
  "VAFB_Arguello"    = "arguello",   "Strauss" = "strauss",
  "VAFB_Lion_Head"   = "lionhead",   "Dangermond" = "Dangermond"
)

# The 13 excluded sample IDs have untrustworthy coordinates.
BAD_COORD_IDS <- c("4548", "4549", "4550", "4551", "4552", "4553", "4554",
                    "4555", "4556", "4557", "4559", "4561", "4562")
field <- read.delim(file.path(BASE_DIR, "field_data.txt"), header = TRUE, sep = "\t")
field$SampleID <- as.character(field$SampleID)
field <- field[!field$SampleID %in% BAD_COORD_IDS, ]
centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field, FUN = mean)
pop_coords <- as.matrix(centroids[, c("Long", "Lat")])
rownames(pop_coords) <- name_map[centroids$Clusters]

# Tranquillon Mtn individuals, same set and order as the OHANA K=6 Q-matrix.
individuals <- read.delim(file.path(WIND_DIR, "strauss_individual", "strauss_individuals.txt"), header = TRUE)
ind_coords <- as.matrix(individuals[, c("lon", "lat")])
rownames(ind_coords) <- individuals$id

# ---------- 3. Pairwise least-cost wind travel times ----------

write_hours <- function(m, ids, file) {
  dimnames(m) <- list(ids, ids)
  write.table(signif(m, 6), file.path(WS_DIR, file), sep = "\t", quote = FALSE)
  cat(sprintf("Saved %d x %d matrix -> %s\n", nrow(m), ncol(m), file))
  m
}

for (clim in c("hrrr", "windninja")) {
  cat(sprintf("\n=== Population centroids, %s wind ===\n", clim))
  hrs <- write_hours(pairwise_least_cost(roses[[clim]], pop_coords),
                     rownames(pop_coords), sprintf("pop_hours_%s.txt", clim))
  print(round(hrs, 1))
}

cat("\n=== Tranquillon Mtn individuals, WindNinja (100 m Tranquillon run) wind ===\n")
t0 <- Sys.time()
hrs <- write_hours(pairwise_least_cost(roses[["windninja_tranquillon"]], ind_coords),
                   rownames(ind_coords), "strauss_hours_windninja.txt")
cat(sprintf("[%.0fs]; travel time range %.2f-%.2f h, NA pairs: %d\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            min(hrs[row(hrs) != col(hrs)], na.rm = TRUE), max(hrs, na.rm = TRUE), sum(is.na(hrs))))

