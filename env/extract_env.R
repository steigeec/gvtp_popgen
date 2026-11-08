#!/usr/bin/env Rscript
# Extract environmental raster values at sample coordinates for env_var.txt.
# Live code adds fog_CV and fog_detrend (2000-2022 MODIS fog) and writes
# env_var_fogCV_fogDetrend_NEW.txt; env_var.txt itself is not modified.
# The commented-out blocks were run once to add the other columns named in
# their headers (fog_germ_yr, BIO*_germ_yr, BIO12/15/18, BIO15_CV/detrend,
# Kling_detrend); uncomment to regenerate them.
# Environment: R library BASE_DIR/software/R_libs
# Usage: Rscript extract_env.R

BASE_DIR <- "/private/groups/meyerlab/esteiger"

.libPaths(file.path(BASE_DIR, "software/R_libs"))
library(terra)
library(ggplot2)
library(patchwork)

setwd(BASE_DIR)

# ==============================================================================
# Fog extraction: mean fog raster per sampling year (fog_germ_yr)
# ==============================================================================
# The fog rasters cover Jun-Sep only, so the Oct (yr-1) to Jun (yr) window
# picks up just the June raster of the sampling year (checked on 60 Tranquillon
# Mtn plants: fog_germ_yr equals the June value exactly). This is intended:
# plants were sampled in May-June, so June is the fog month overlapping the
# sampled cohort's season up to sampling.

# ---- Build per-year mean fog rasters (Oct yr-1 through Jun yr) ---------------
# fog_dir   <- "env/fog"
# fog_files <- list.files(fog_dir, pattern="^\\d{4}-\\d{2}\\.tif$", full.names=TRUE)
# years     <- c(2020, 2021, 2022)

# Parse dates from fog filenames (expected format: YYYY-MM.tif)
# fog_file_dates <- as.Date(paste0(
#   sub("\\.tif$", "", basename(fog_files)), "-01"))

# fog_mean <- list()
# for (yr in years) {
#   yr_chr  <- as.character(yr)
#   prev_yr <- yr - 1

#   # Keep only fog files within Oct (yr-1) through Jun (yr)
#   in_range <- (format(fog_file_dates, "%Y") == as.character(prev_yr) &
#                  as.integer(format(fog_file_dates, "%m")) >= 10) |
#               (format(fog_file_dates, "%Y") == yr_chr &
#                  as.integer(format(fog_file_dates, "%m")) <= 6)

#   yr_files <- fog_files[in_range]
#   if (length(yr_files) == 0) {
#     warning(sprintf("No fog files found for germ year %d (Oct %d – Jun %d)",
#                     yr, prev_yr, yr))
#     next
#   }
#   stack <- rast(yr_files)
#   fog_mean[[yr_chr]] <- mean(stack, na.rm=TRUE)
#   cat(sprintf("Germ year %d (Oct %d – Jun %d): averaged %d fog rasters (%s)\n",
#               yr, prev_yr, yr, length(yr_files),
#               paste(basename(yr_files), collapse=", ")))
# }

# ---- Load sample coordinates -------------------------------------------------
env <- read.table("env_var.txt", header=TRUE, sep="\t")
cat(sprintf("\nenv_var.txt: %d rows, %d columns\n", nrow(env), ncol(env)))

valid <- !is.na(env$Lat) & !is.na(env$Long) &
         env$Lat != -999999 & env$Long != -999999
cat(sprintf("Samples with valid coordinates: %d / %d\n", sum(valid), nrow(env)))

# ---- Extract fog for each sample's sampling year ----------------------------
# env$fog_germ_yr <- NA_real_

# for (yr in years) {
#   yr_chr  <- as.character(yr)
#   if (is.null(fog_mean[[yr_chr]])) next

#   idx <- valid & env$Year == yr
#   if (sum(idx) == 0) next

#   pts <- vect(env[idx, ], geom=c("Long","Lat"), crs="EPSG:4326")
#   ext_vals <- extract(fog_mean[[yr_chr]], pts)
#   env$fog_germ_yr[idx] <- ext_vals[, 2]

#   cat(sprintf("Year %d: extracted fog for %d samples (%d NA)\n",
#               yr, sum(idx), sum(is.na(ext_vals[, 2]))))
# }

# cat(sprintf("\nfog_germ_yr summary:\n"))
# print(summary(env$fog_germ_yr))

# ---- Write updated env_var.txt -----------------------------------------------
# write.table(env, "env_var.txt", sep="\t", quote=FALSE, row.names=FALSE)
# cat("\nUpdated env_var.txt written with fog_germ_yr column.\n")

# ==============================================================================
# BIO12_germ_yr, BIO15_germ_yr, BIO18_germ_yr
# Source: PRISM monthly 800m precipitation
# Water year = Oct yr-1 through Sep yr (12 months)
#   BIO12: total water-year precipitation (mm)
#   BIO15: precipitation seasonality, CV x 100 of the 12 monthly values
#   BIO18: precipitation of warmest quarter (Jun + Jul + Aug of sampling year)
# ==============================================================================

# prism_dir <- "env/PRISM"

# Helper: locate 800m PRISM .tif for a given YYYYMM
# find_prism_tif <- function(yyyymm) {
#   files <- list.files(file.path(prism_dir, yyyymm),
#                       pattern = "_30s_.*\\.tif$", full.names = TRUE)
#   files <- files[!grepl("\\.aux\\.xml$", files)]
#   if (length(files) == 0) return(NULL)
#   files[1]
# }

# ---- Extract PRISM water-year variables per sampling year --------------------
# env$BIO12_germ_yr <- NA_real_
# env$BIO15_germ_yr <- NA_real_
# env$BIO18_germ_yr <- NA_real_

# for (yr in years) {
#   idx <- valid & env$Year == yr
#   if (sum(idx) == 0) next

#   wy_yyyymm    <- c(sprintf("%04d%02d", yr - 1, 10:12),
#                     sprintf("%04d%02d", yr,      1:9))
#   bio18_yyyymm <- sprintf("%04d%02d", yr, 6:8)   # Jun–Aug of sampling year

#   tif_paths <- sapply(wy_yyyymm, find_prism_tif)
#   missing   <- wy_yyyymm[sapply(tif_paths, is.null)]
#   if (length(missing) > 0) {
#     warning(sprintf("Missing PRISM raster(s) for water year %d: %s",
#                     yr, paste(missing, collapse = ", ")))
#     next
#   }

#   stack        <- rast(unlist(tif_paths))
#   names(stack) <- wy_yyyymm
#   pts_yr       <- vect(env[idx, ], geom = c("Long", "Lat"), crs = "EPSG:4326")
#   pts_prj      <- project(pts_yr, crs(stack))
#   vals         <- as.matrix(extract(stack, pts_prj)[, -1])

#   env$BIO12_germ_yr[idx] <- rowSums(vals)
#   env$BIO15_germ_yr[idx] <- apply(vals, 1, function(x) {
#     if (any(is.na(x))) return(NA_real_)
#     (sd(x) / mean(x)) * 100
#   })
#   env$BIO18_germ_yr[idx] <- rowSums(vals[, bio18_yyyymm, drop = FALSE])

#   cat(sprintf("Water year %d: extracted PRISM variables for %d samples (%d NA in BIO12)\n",
#               yr, sum(idx), sum(is.na(env$BIO12_germ_yr[idx]))))
# }

# cat("\nBIO12_germ_yr summary:\n")
# print(summary(env$BIO12_germ_yr))
# cat("\nBIO15_germ_yr summary:\n")
# print(summary(env$BIO15_germ_yr))
# cat("\nBIO18_germ_yr summary:\n")
# print(summary(env$BIO18_germ_yr))

# ---- Write updated env_var.txt -----------------------------------------------
# write.table(env, "env_var.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# cat("\nUpdated env_var.txt written with BIO12_germ_yr, BIO15_germ_yr, BIO18_germ_yr columns.\n")

# ==============================================================================
# BIO12, BIO15, BIO18: PRISM averages over water years 2000-2022 (the fog period)
#   BIO12 : mean annual precipitation (mm/yr)
#   BIO15 : precipitation seasonality, CV x 100 of the 12 long-term monthly means
#   BIO18 : mean warm-season precipitation (Jun + Jul + Aug, mm)
# ==============================================================================

# bio_wy_range <- 2000:2022
# n_wy         <- length(bio_wy_range)

# cat("\n--- Computing PRISM-based BIO12/BIO15/BIO18 (WY2000-2022 averages) ---\n")

# pts_bio <- vect(env[valid, ], geom = c("Long", "Lat"), crs = "EPSG:4326")
# n_pts   <- sum(valid)

# Accumulators: rows = valid samples, cols = water years / month positions
# bio12_mat  <- matrix(NA_real_, nrow = n_pts, ncol = n_wy)
# bio18_mat  <- matrix(NA_real_, nrow = n_pts, ncol = n_wy)
# month_sums <- matrix(0,        nrow = n_pts, ncol = 12)   # sum over years per WY month position
# n_wy_done  <- 0

# for (i in seq_along(bio_wy_range)) {
#   yr           <- bio_wy_range[i]
#   wy_yyyymm    <- c(sprintf("%04d%02d", yr - 1, 10:12),
#                     sprintf("%04d%02d", yr,      1:9))
#   bio18_yyyymm <- sprintf("%04d%02d", yr, 6:8)   # Jun–Aug of yr

#   tif_paths <- sapply(wy_yyyymm, find_prism_tif)
#   missing   <- wy_yyyymm[sapply(tif_paths, is.null)]
#   if (length(missing) > 0) {
#     warning(sprintf("Skipping WY%d — missing PRISM raster(s): %s",
#                     yr, paste(missing, collapse = ", ")))
#     next
#   }

#   stack        <- rast(unlist(tif_paths))
#   names(stack) <- wy_yyyymm
#   pts_prj      <- project(pts_bio, crs(stack))
#   vals         <- as.matrix(extract(stack, pts_prj)[, -1])

#   bio12_mat[, i] <- rowSums(vals)
#   bio18_mat[, i] <- rowSums(vals[, bio18_yyyymm, drop = FALSE])
#   month_sums     <- month_sums + vals   # accumulate monthly totals
#   n_wy_done      <- n_wy_done + 1
#   cat(sprintf("  WY%d done\n", yr))
# }

# Long-term monthly means (one per WY month position across completed years)
# month_ltm <- month_sums / n_wy_done

# env$BIO12[valid] <- rowMeans(bio12_mat, na.rm = TRUE)
# env$BIO15[valid] <- apply(month_ltm, 1, function(x) {
#   if (any(is.na(x))) return(NA_real_)
#   (sd(x) / mean(x)) * 100
# })
# env$BIO18[valid] <- rowMeans(bio18_mat, na.rm = TRUE)

# cat(sprintf("\nWater years completed: %d / %d\n", n_wy_done, n_wy))
# cat(sprintf("BIO12 (PRISM WY2000-2022 mean): %.1f ± %.1f mm\n",
#             mean(env$BIO12[valid], na.rm = TRUE), sd(env$BIO12[valid], na.rm = TRUE)))
# cat(sprintf("BIO15 (PRISM WY2000-2022):      %.1f ± %.1f\n",
#             mean(env$BIO15[valid], na.rm = TRUE), sd(env$BIO15[valid], na.rm = TRUE)))
# cat(sprintf("BIO18 (PRISM WY2000-2022 mean): %.1f ± %.1f mm\n",
#             mean(env$BIO18[valid], na.rm = TRUE), sd(env$BIO18[valid], na.rm = TRUE)))

# write.table(env, "env_var.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# cat("Updated env_var.txt written with PRISM-based BIO12, BIO15, BIO18.\n")

# ==============================================================================
# FOG VARIABILITY: fog_CV and fog_detrend, 2000-2022
# Annual Jun-Sep mean fog per point; its SD is the existing fog_SD column.
#   fog_CV      : SD / mean x 100 across the 23 annual values; normalizing by
#                 each site's mean fog (as for BIO15) avoids confounding with mean fog
#   fog_detrend : SD of residuals after removing a linear year trend, so a
#                 steady long-term change does not inflate the variability measure
# ==============================================================================

fog_sd_dir   <- "env/fog"
fog_sd_files <- list.files(fog_sd_dir, pattern="^\\d{4}-\\d{2}\\.tif$", full.names=TRUE)
fog_sd_dates <- as.Date(paste0(
  sub("\\.tif$", "", basename(fog_sd_files)), "-01"))

fog_sd_years <- 2000:2022

pts_fogsd      <- vect(env[valid, ], geom=c("Long","Lat"), crs="EPSG:4326")
fog_annual_mat <- matrix(NA_real_, nrow=sum(valid), ncol=length(fog_sd_years))

for (i in seq_along(fog_sd_years)) {
  yr <- fog_sd_years[i]
  yr_files <- fog_sd_files[format(fog_sd_dates, "%Y") == as.character(yr)]
  if (length(yr_files) == 0) {
    warning(sprintf("No fog files found for year %d", yr))
    next
  }
  # extract each month separately: not all monthly rasters share a grid (e.g. 2022-09.tif)
  month_vals <- sapply(yr_files, function(f) extract(rast(f), pts_fogsd)[, 2])
  fog_annual_mat[, i] <- rowMeans(month_vals, na.rm=TRUE)
  cat(sprintf("Fog year %d: averaged %d rasters (%s)\n", yr, length(yr_files),
              paste(basename(yr_files), collapse=", ")))
}

# ---- Check against the existing fog_SD column ------------------------------
fog_sd_recomputed <- apply(fog_annual_mat, 1, sd, na.rm=TRUE)
existing_fog_sd    <- env$fog_SD[valid]
sd_diff <- abs(fog_sd_recomputed - existing_fog_sd)
cat(sprintf(
  "\nfog_SD verification: max abs diff vs existing column = %.10f (%d / %d within 1e-8)\n",
  max(sd_diff, na.rm = TRUE), sum(sd_diff < 1e-8, na.rm = TRUE), sum(valid)))

# ---- fog_CV: SD / mean x 100 ------------------------------------------------
fog_mean_annual <- rowMeans(fog_annual_mat, na.rm = TRUE)
fog_cv <- (fog_sd_recomputed / fog_mean_annual) * 100

# ---- fog_detrend: SD of residuals after removing a linear year trend -------
# Uses annual means rather than regressing raw monthly values on year +
# factor(month): the monthly version shifted the population ranking (Spearman
# rho = 0.7 across the 5 TajD populations; Gaviota went from lowest to
# 3rd-lowest variability) enough to change which population anchors the
# fog-variability mixed model.
fog_detrend <- apply(fog_annual_mat, 1, function(x) {
  ok <- !is.na(x)
  if (sum(ok) < 3) return(NA_real_)
  sd(residuals(lm(x[ok] ~ fog_sd_years[ok])))
})

env$fog_CV      <- NA_real_
env$fog_detrend <- NA_real_
env$fog_CV[valid]      <- fog_cv
env$fog_detrend[valid] <- fog_detrend

cat("\nfog_CV summary (SD/mean x 100 across 2000-2022 of Jun-Sep mean fog):\n")
print(summary(env$fog_CV))
cat("\nfog_detrend summary (SD of residuals after removing a linear year trend):\n")
print(summary(env$fog_detrend))

# ---- Write to a new file; env_var.txt is not modified ----------------------
out_file <- "env_var_fogCV_fogDetrend_NEW.txt"
write.table(env, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nWritten %s (added fog_CV, fog_detrend; nothing else modified).\n", out_file))

# # ==============================================================================
# # BIO15 VARIABILITY: BIO15_CV and BIO15_detrend, full PRISM record
# # BIO15 computed per water year over the full gap-free PRISM record
# # (1927-2024, not capped to the 2000-2022 fog window), then per point:
# #   BIO15_CV      : SD / mean x 100 across the annual BIO15 values
# #   BIO15_detrend : SD of residuals after removing a linear year trend
# # ==============================================================================

# prism_dir <- "env/PRISM"

# # Helper: locate 800m PRISM .tif for a given YYYYMM
# find_prism_tif <- function(yyyymm) {
#   files <- list.files(file.path(prism_dir, yyyymm),
#                       pattern = "_30s_.*\\.tif$", full.names = TRUE)
#   files <- files[!grepl("\\.aux\\.xml$", files)]
#   if (length(files) == 0) return(NULL)
#   files[1]
# }

# bio15_years <- 1927:2024

# pts_bio15        <- vect(env[valid, ], geom = c("Long", "Lat"), crs = "EPSG:4326")
# bio15_annual_mat <- matrix(NA_real_, nrow = sum(valid), ncol = length(bio15_years))

# for (i in seq_along(bio15_years)) {
#   yr        <- bio15_years[i]
#   wy_yyyymm <- c(sprintf("%04d%02d", yr - 1, 10:12),
#                  sprintf("%04d%02d", yr,      1:9))

#   tif_paths <- sapply(wy_yyyymm, find_prism_tif)
#   missing   <- wy_yyyymm[sapply(tif_paths, is.null)]
#   if (length(missing) > 0) {
#     warning(sprintf("Skipping WY%d for BIO15 — missing PRISM raster(s): %s",
#                     yr, paste(missing, collapse = ", ")))
#     next
#   }

#   # some PRISM tifs are truncated on disk (e.g. 194811); skip that water year
#   vals <- tryCatch({
#     stack   <- rast(unlist(tif_paths))
#     pts_prj <- project(pts_bio15, crs(stack))
#     as.matrix(extract(stack, pts_prj)[, -1])
#   }, error = function(e) {
#     warning(sprintf("Skipping WY%d for BIO15 — error reading PRISM raster(s): %s",
#                     yr, conditionMessage(e)))
#     NULL
#   })
#   if (is.null(vals)) next

#   bio15_annual_mat[, i] <- apply(vals, 1, function(x) {
#     if (any(is.na(x))) return(NA_real_)
#     (sd(x) / mean(x)) * 100
#   })
#   cat(sprintf("BIO15 water year %d: computed annual BIO15 for %d samples\n",
#               yr, sum(valid)))
# }

# bio15_mean_annual <- rowMeans(bio15_annual_mat, na.rm = TRUE)
# bio15_sd_annual   <- apply(bio15_annual_mat, 1, sd, na.rm = TRUE)
# bio15_cv          <- (bio15_sd_annual / bio15_mean_annual) * 100

# bio15_detrend <- apply(bio15_annual_mat, 1, function(x) {
#   ok <- !is.na(x)
#   if (sum(ok) < 3) return(NA_real_)
#   sd(residuals(lm(x[ok] ~ bio15_years[ok])))
# })

# env$BIO15_CV      <- NA_real_
# env$BIO15_detrend <- NA_real_
# env$BIO15_CV[valid]      <- bio15_cv
# env$BIO15_detrend[valid] <- bio15_detrend

# cat(sprintf("\nBIO15_CV summary (SD/mean x 100 across %d-%d of annual BIO15):\n",
#             min(bio15_years), max(bio15_years)))
# print(summary(env$BIO15_CV))
# cat("\nBIO15_detrend summary (SD of residuals after removing a linear year trend):\n")
# print(summary(env$BIO15_detrend))

# out_file2 <- "env_var_BIO15CV_BIO15detrend_NEW.txt"
# write.table(env, out_file2, sep = "\t", quote = FALSE, row.names = FALSE)
# cat(sprintf("\nWritten %s (added BIO15_CV, BIO15_detrend; nothing else modified).\n", out_file2))

# # ==============================================================================
# # KLING FOG VARIABILITY: Kling_detrend, full long-term record (1950-2014)
# # env/fog/Kling/timeseries.tif has one band per year.
# #   Kling_detrend : SD of residuals after removing a linear year trend
# # ==============================================================================

# kling_r     <- rast("env/fog/Kling/timeseries.tif")
# kling_years <- as.integer(names(kling_r))

# pts_kling      <- vect(env[valid, ], geom = c("Long", "Lat"), crs = "EPSG:4326")
# pts_kling_prj  <- project(pts_kling, crs(kling_r))
# kling_annual_mat <- as.matrix(extract(kling_r, pts_kling_prj)[, -1])

# kling_detrend <- apply(kling_annual_mat, 1, function(x) {
#   ok <- !is.na(x)
#   if (sum(ok) < 3) return(NA_real_)
#   sd(residuals(lm(x[ok] ~ kling_years[ok])))
# })

# env$Kling_detrend <- NA_real_
# env$Kling_detrend[valid] <- kling_detrend

# cat(sprintf("\nKling_detrend summary (SD of residuals after removing a linear year trend, %d-%d):\n",
#             min(kling_years), max(kling_years)))
# print(summary(env$Kling_detrend))

# out_file3 <- "env_var_Kling_detrend_NEW.txt"
# write.table(env, out_file3, sep = "\t", quote = FALSE, row.names = FALSE)
# cat(sprintf("\nWritten %s (added Kling_detrend; nothing else modified).\n", out_file3))
