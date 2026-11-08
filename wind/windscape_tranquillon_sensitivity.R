## Tranquillon Mtn isolation-by-wind tests for each WindNinja sensitivity
## variant (run_windninja_tranquillon_sensitivity.sh) and the full
## configuration (run_windninja_landfire_batch.sh tranquillon), using the same
## windscape steps and tests as windscape_tests.R.
##
## Output: BASE_DIR/wind/windscape/sensitivity/<variant>_strauss_hours.txt,
##   <variant>_results.txt, tranquillon_sensitivity_results.txt (combined)
## R libraries: windscape 2.0.0 and terra 1.9-50 (CRAN) in BASE_DIR/software/R_libs_windscape,
##   installed with install.packages(c("windscape", "terra"), lib = <that dir>);
##   geosphere from BASE_DIR/software/R_libs
## Usage: Rscript windscape_tranquillon_sensitivity.R base mesh100 landfire diurnal stability full
##   Variants can run as separate parallel calls; each call rebuilds the combined table.

BASE_DIR <- "/private/groups/meyerlab/esteiger"
.libPaths(c(file.path(BASE_DIR, "software/R_libs_windscape"),
            file.path(BASE_DIR, "software/R_libs/")))
library(windscape)
library(terra)
library(geosphere)

WIND_DIR <- file.path(BASE_DIR, "wind")
WS_DIR <- file.path(WIND_DIR, "windscape")
SENS_OUT <- file.path(WS_DIR, "sensitivity")
dir.create(SENS_OUT, showWarnings = FALSE)
NPERM <- 99999
TRANQ_EXT <- ext(-120.558, -120.496, 34.559, 34.609)
TRANQ_RES <- 0.001

variant_dir <- function(v) {
  if (v == "full") file.path(WIND_DIR, "windninja_landfire", "tranquillon")
  else file.path(WIND_DIR, "windninja_landfire", "sensitivity", v)
}

# UTC timestamps of the bloom-season sample, from the HRRR series' layer names.
stamp <- sub("^u_", "", grep("^u_", names(rast(file.path(WS_DIR, "hrrr_wind_series.tif"))), value = TRUE))

# ---------- Shared genetic / geographic / environmental inputs ----------

individuals <- read.delim(file.path(WIND_DIR, "strauss_individual", "strauss_individuals.txt"), header = TRUE)
ind_coords <- as.matrix(individuals[, c("lon", "lat")])
rownames(ind_coords) <- individuals$id

bam_ids <- sub("_rmdup\\.bam$", "", readLines(file.path(BASE_DIR, "mapping/merged/bamlist_strauss.txt")))
fst_all <- read.delim(file.path(BASE_DIR, "structure/Fst/individuals/pairwise_Fst.txt"), header = TRUE)
fst_all$sample1 <- as.character(fst_all$sample1); fst_all$sample2 <- as.character(fst_all$sample2)
fst_s <- fst_all[fst_all$sample1 %in% bam_ids & fst_all$sample2 %in% bam_ids, ]
fst_s <- fst_s[fst_s$sample1 != "6432" & fst_s$sample2 != "6432", ]  # duplicate-sampled plant
ids <- sort(unique(c(fst_s$sample1, fst_s$sample2)))
fst_mat <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
fst_mat[cbind(fst_s$sample1, fst_s$sample2)] <- fst_s$Fst_weighted
fst_mat[cbind(fst_s$sample2, fst_s$sample1)] <- fst_s$Fst_weighted
fst_lin <- fst_mat / (1 - fst_mat)

field <- read.delim(file.path(BASE_DIR, "field_data.txt"), header = TRUE, sep = "\t")
field$SampleID <- as.character(field$SampleID)
geo_m <- distm(as.matrix(field[match(ids, field$SampleID), c("Long", "Lat")]), fun = distHaversine)
geo_m[geo_m == 0 & row(geo_m) != col(geo_m)] <- 1  # one coincident pair, floored at 1 m
geo_log <- log(geo_m / 1000)
diag(geo_log) <- 0
env_raw <- read.delim(file.path(BASE_DIR, "env_var.txt"), header = TRUE)
env_raw$SampleID <- as.character(env_raw$SampleID)
env_dist <- as.matrix(dist(scale(env_raw[match(ids, env_raw$SampleID), c("elevation", "BIO18", "fog")])))

# ---------- Per-variant pipeline ----------

run_variant <- function(v) {
  t0 <- Sys.time()
  vel_files <- list.files(variant_dir(v), pattern = "_vel\\.asc$", full.names = TRUE)
  m <- regmatches(basename(vel_files), regexec("_(\\d{2})-(\\d{2})-(\\d{4})_(\\d{2})(\\d{2})_\\d+m.*_vel\\.asc$", basename(vel_files)))
  local <- as.POSIXct(sapply(m, function(x) sprintf("%s-%s-%s %s:%s", x[4], x[2], x[3], x[5], x[6])),
                      format = "%Y-%m-%d %H:%M", tz = "America/Los_Angeles")
  vel_by_stamp <- setNames(vel_files, format(local, "%Y%m%dT%H", tz = "UTC"))
  stopifnot(all(stamp %in% names(vel_by_stamp)))

  tmpl <- rast(TRANQ_EXT, resolution = TRANQ_RES, crs = "EPSG:4326")
  u <- v_ <- vector("list", length(stamp))
  for (k in seq_along(stamp)) {
    vel <- rast(vel_by_stamp[stamp[k]])
    rad <- rast(sub("_vel\\.asc$", "_ang.asc", vel_by_stamp[stamp[k]])) * pi / 180
    uv <- c(-vel * sin(rad), -vel * cos(rad))
    if (is.na(crs(uv)) || crs(uv) == "") crs(uv) <- "EPSG:32610"
    uv <- project(uv, tmpl, method = "bilinear")
    u[[k]] <- uv[[1]]; v_[[k]] <- uv[[2]]
  }
  series <- c(rast(u), rast(v_))
  stopifnot(!anyNA(values(series[[1]])))
  rose <- wind_rose(wind_series(series), trans = 1)
  hours <- pairwise_least_cost(rose, ind_coords)
  dimnames(hours) <- list(rownames(ind_coords), rownames(ind_coords))
  write.table(signif(hours, 6), file.path(SENS_OUT, sprintf("%s_strauss_hours.txt", v)), sep = "\t", quote = FALSE)

  h <- hours[ids, ids]
  off <- row(hours) != col(hours)
  h[h == 0 & row(h) != col(h)] <- min(hours[off & hours > 0])  # coincident individuals
  wind_dist <- -log(pairwise_means(1 / h))
  diag(wind_dist) <- 0

  set.seed(42)
  geo <- mantel_test(fst_lin, wind_dist, list(geo_log), nperm = NPERM)
  env <- mantel_test(fst_lin, wind_dist, list(env_dist), nperm = NPERM)
  both <- mantel_test(fst_lin, wind_dist, list(geo_log, env_dist), nperm = NPERM)
  cat(sprintf("%-10s | geo r = %.4f p = %.4f | env r = %.4f p = %.4f | geo+env r = %.4f p = %.4f [%.0fs]\n",
              v, geo$stat, geo$p.value, env$stat, env$p.value, both$stat, both$p.value,
              as.numeric(Sys.time() - t0, units = "secs")))
  data.frame(variant = v, n_timesteps = length(stamp),
             r_geo = geo$stat, p_geo = geo$p.value, r_env = env$stat, p_env = env$p.value,
             r_geo_env = both$stat, p_geo_env = both$p.value)
}

for (v in commandArgs(trailingOnly = TRUE)) {
  write.table(run_variant(v), file.path(SENS_OUT, sprintf("%s_results.txt", v)), sep = "\t", quote = FALSE, row.names = FALSE)
}
order <- c("base", "mesh100", "landfire", "diurnal", "stability", "full")
files <- file.path(SENS_OUT, sprintf("%s_results.txt", order))
res <- do.call(rbind, lapply(files[file.exists(files)], read.delim))
write.table(res, file.path(SENS_OUT, "tranquillon_sensitivity_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, digits = 4)
