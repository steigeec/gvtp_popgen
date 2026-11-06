# extract_pheno_env_data.R
# Extract environmental data for phenology samples relative to sampling year.
# Output: emma_gee_results_new.csv, env_wi_range.csv,
#         env_selection_correlation_heatmap_25yr_modeling_vars.jpg (Section 13)
#
# For each sample, up to three time-averaged variants are computed,
# depending on what the underlying data source supports:
#   _germyr : value for the specific year the plant was collected ("Year")
#   _25yr   : mean over the fixed reference window common to (most) sources
#             — 2000-2025 for fog (calendar years), 2000-2024 water years
#             for BIO12/15/18 (PRISM archive ends Sep 2024)
#   _<N>yr  : mean over the long-term window, sized and dated to match
#             Kling fog's real observed range (1950-2014, N=65) exactly —
#             BIO12/15/18 are averaged over those same water years, not
#             the full 1927-2024 PRISM archive, so the "long" fog and
#             bioclim columns are directly comparable over identical years.
#   CHILI is a static topographic layer (no time dimension) and gets a
#   single column with no suffix.
#   Kling has no "_25yr" column: its data (1950-2014) doesn't reach 2025,
#   so the fixed reference window isn't fully covered.
#   fog has no "_<N>yr" column: its full available range (2000-2025) is
#   identical to "_25yr", so a third variant would be redundant.
#
# Variables:
#   fog:   Mean monthly FLCC (fog/low cloud cover days) June-September
#           Source: MODIS Monthly FLCC Rasters 2000-2025
#           (Werner et al. 2022, doi:10.25675/10217/235754;
#            2023-2025 rasters provided by Zackary Werner)
#   CHILI: Continuous Heat-Insolation Load Index (static topographic layer)
#   BIO12: Total annual precipitation over the water year (Oct yr-1 – Sep yr)
#   BIO15: Precipitation seasonality — CV of 12 monthly values
#   BIO18: Precipitation of the warmest quarter (Jun + Jul + Aug of the year)
#           Source: PRISM monthly 800m precipitation
#   Kling: Annual (Jun-Sep mean) fog frequency from the Kling lab's observed
#           time series (real per-year data — supersedes the earlier
#           trend+normal reconstruction)
#           Source: Kling/timeseries.tif, 65 annual layers, 1950-2014
#   CWD:   Climatic water deficit (mm/water-year), summed monthly.
#           Only a "_25yr" column (see Section 12) — the BCM archive on
#           disk covers WY2000-2022 (not 2023/2024), so CWD_25yr is
#           actually a 23-water-year mean, the closest available analog
#           to the other "_25yr" columns.
#           Source: BCMv8 statewide monthly CWD, EPSG:3310, 270 m
#           (see env/extract_BCM_CWD.R for the ScienceBase source URLs)

library(terra)
library(dplyr)
library(ggplot2)

# ============================================================
# Paths and parameters
# ============================================================

coords_file   <- "/private/groups/meyerlab/esteiger/pheno/all_coord_gvtp.csv"
fog_dir       <- "/private/groups/meyerlab/esteiger/env/fog/"
chili_file    <- "/private/groups/meyerlab/esteiger/env/US_CHILI.tif"
prism_dir     <- "/private/groups/meyerlab/esteiger/env/PRISM/"
kling_ts_file <- "/private/groups/meyerlab/esteiger/env/fog/Kling/timeseries.tif"
kling_ts_zip  <- "/private/groups/meyerlab/esteiger/env/fog/Kling/timeseries.tif.zip"
output_file   <- "/private/groups/meyerlab/esteiger/pheno/emma_gee_results_final.csv"
plots_dir     <- "/private/groups/meyerlab/esteiger/pheno/"
pheno_file    <- "/private/groups/meyerlab/esteiger/pheno/gvtp_phenotype_13feb2026.tsv"

# Crop extent for PRISM rasters: widened 2026-10-27 to California-wide
# (was ext(-120.75,-120.1,34.4,35.2), the tiny Pt Conception study bbox) so
# the archive can also serve the standalone CA summer-precip trend map
# (env/download_crop_prism_ca.R) without a second copy of the data. Still a
# strict superset of every sample coordinate, so nothing here changes for
# the point-extraction logic below -- this just keeps Section 3b's
# crop_ppt_archive_inplace() (only runs if RUN_FULL_EXTRACTION <- TRUE) from
# cropping the CA-wide files back down to the old tiny bbox if ever re-run.
bbox_extent   <- ext(-124.6, -113.9, 32.4, 42.2)

fog_months    <- 6:9   # June through September
fog_year_min  <- 2000
fog_year_max  <- 2025

# Fixed reference window shared across variables ("_25yr" columns)
ref_fog_years <- fog_year_min:fog_year_max   # 2000-2025 (26 calendar years)

# Full historical PRISM archive: Oct 1926 – Sep 2025 -> water years 1927-2025
# (WY2025 added 2026-09-07 specifically to cover 2025/2026-dated iNat
# observations, whose germination years exceeded the archive's prior
# WY2024 ceiling and were excluded from BIO18_modern/BIO15_modern entirely)
prism_wy_first <- 1927
prism_wy_last  <- 2025
ref_bio_years  <- 2000:prism_wy_last         # 2000-2025 (26 water years)

# Full dataset zip (MODIS Monthly FLCC Rasters 2000-2022, CSU Libraries)
# Note: 2023-2025 rasters were provided by Zackary Werner and must be present
# in fog_dir before running.
fog_zip_url <- "https://mountainscholar.org/bitstreams/9c4dea5c-bf1c-4d30-bd9e-37bbab77d855/download"

# ============================================================
# Sections 1-11 (fog/PRISM/CHILI/Kling download, extraction, bioclim
# analogs, collinearity/PCA, and the env_wi_range.csv full-year-range
# export) already ran to completion — their outputs are cached on disk
# (emma_gee_results_final.csv, env_wi_range.csv, the diagnostic jpgs).
# Disabled wholesale (rather than re-run) so this script now only adds
# the new BCM CWD_25yr column — see "Section 12" at the bottom.
# Flip to TRUE and re-run from scratch if the underlying rasters change.
# ============================================================

RUN_FULL_EXTRACTION <- FALSE

if (RUN_FULL_EXTRACTION) {

# ============================================================
# 1. Load coordinates
# ============================================================

coords <- read.csv(coords_file)
cat("Loaded", nrow(coords), "samples\n")

sample_years <- sort(unique(coords$Year))
cat("Sample years:", paste(sample_years, collapse = ", "), "\n")

years_out_range <- sample_years[sample_years < fog_year_min | sample_years > fog_year_max]
if (length(years_out_range) > 0) {
  warning(
    "The following sample years fall outside the fog dataset range (",
    fog_year_min, "-", fog_year_max,
    ") and will have NA fog_germyr values: ",
    paste(years_out_range, collapse = ", ")
  )
}

# ============================================================
# 2. Check for missing fog files; download and extract if needed
#    (Full 2000-2025 range required for the fixed "_25yr" average,
#    regardless of which years the samples themselves fall in.)
# ============================================================

needed <- expand.grid(year = ref_fog_years, month = fog_months)
needed$filename <- sprintf("%d-%02d.tif", needed$year, needed$month)
needed$path     <- file.path(fog_dir, needed$filename)
needed$exists   <- file.exists(needed$path)

missing <- needed[!needed$exists, ]

if (nrow(missing) == 0) {
  cat("All required fog files are already present.\n\n")
} else {
  cat("\nMissing fog files (", nrow(missing), "):\n", sep = "")
  cat(paste(" ", missing$filename, collapse = "\n"), "\n\n")

  zip_path <- tempfile(fileext = ".zip")
  cat("Downloading MODIS Monthly FLCC Rasters dataset (~87 MB)...\n")
  download.file(fog_zip_url, zip_path, mode = "wb", quiet = FALSE)

  cat("Listing zip contents...\n")
  zip_contents <- unzip(zip_path, list = TRUE)$Name

  cat("Extracting missing files...\n")
  for (f in missing$filename) {
    match <- zip_contents[grepl(f, zip_contents, fixed = TRUE)]
    if (length(match) == 0) {
      warning("File '", f, "' not found in zip archive — skipping.")
      next
    }
    unzip(zip_path, files = match[1], exdir = fog_dir, junkpaths = TRUE)
    cat("  Extracted:", f, "\n")
  }

  unlink(zip_path)
  cat("Download and extraction complete.\n\n")
}

# ============================================================
# 3. Check for missing PRISM files; download if needed
#    (Full water-year range 1927-2024 required for the "_25yr" and
#    "_98yr" averages.)
# ============================================================

dir.create(prism_dir, showWarnings = FALSE)

find_prism_bil <- function(yyyymm) {
  dir_path <- file.path(prism_dir, yyyymm)
  files    <- list.files(dir_path, pattern = "_30s_.*\\.(tif|bil)$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  files[1]
}

wy_months_needed <- do.call(rbind, lapply(prism_wy_first:prism_wy_last, function(yr) {
  rbind(
    data.frame(year = yr - 1, month = 10:12),
    data.frame(year = yr,     month = 1:9)
  )
}))
wy_months_needed$yyyymm <- sprintf("%04d%02d", wy_months_needed$year, wy_months_needed$month)
unique_yyyymm <- sort(unique(wy_months_needed$yyyymm))

cat("Checking PRISM monthly precipitation files...\n")
missing_prism <- unique_yyyymm[sapply(unique_yyyymm, function(m) is.null(find_prism_bil(m)))]

if (length(missing_prism) == 0) {
  cat("All required PRISM files are already present.\n\n")
} else {
  cat("Downloading", length(missing_prism), "missing PRISM monthly file(s)...\n")
  for (yyyymm in missing_prism) {
    url      <- paste0("https://services.nacse.org/prism/data/get/us/800m/ppt/", yyyymm)
    zip_path <- tempfile(fileext = ".zip")
    cat("  Downloading:", yyyymm, "\n")
    tryCatch({
      download.file(url, zip_path, mode = "wb", quiet = TRUE)
      out_dir <- file.path(prism_dir, yyyymm)
      dir.create(out_dir, showWarnings = FALSE)
      unzip(zip_path, exdir = out_dir)
      unlink(zip_path)
    }, error = function(e) {
      warning("Failed to download PRISM data for ", yyyymm, ": ", e$message)
      unlink(zip_path)
    })
  }
  cat("PRISM download complete.\n\n")
}

# ============================================================
# 3b. Crop existing full-CONUS PRISM ppt archive in place to bbox_extent
#     The archive was originally downloaded at full CONUS extent
#     (~60 MB/file, ~66 GB total for the 1927-2024 range). Every sample
#     coordinate falls within bbox_extent, so crop each file in place to
#     free the unused extent. Idempotent: skips files already cropped.
# ============================================================

crop_ppt_archive_inplace <- function() {
  ppt_files <- list.files(prism_dir, pattern = "^prism_ppt_us_30s_.*\\.tif$",
                           recursive = TRUE, full.names = TRUE)
  if (length(ppt_files) == 0) {
    cat("No PRISM ppt files found to check for cropping.\n\n")
    return(invisible(NULL))
  }
  cat("Checking", length(ppt_files), "PRISM ppt file(s) for cropping to bbox_extent...\n")
  bbox_v <- as.vector(bbox_extent)
  n_cropped <- 0
  for (f in ppt_files) {
    r   <- rast(f)
    e_v <- as.vector(ext(r))
    already_cropped <- e_v["xmin"] >= bbox_v["xmin"] - 1e-6 && e_v["xmax"] <= bbox_v["xmax"] + 1e-6 &&
                        e_v["ymin"] >= bbox_v["ymin"] - 1e-6 && e_v["ymax"] <= bbox_v["ymax"] + 1e-6
    if (already_cropped) next
    r_crop <- crop(r, bbox_extent)
    tmp    <- tempfile(fileext = ".tif")
    writeRaster(r_crop, tmp, overwrite = TRUE)
    rm(r, r_crop)
    file.copy(tmp, f, overwrite = TRUE)
    file.remove(tmp)
    n_cropped <- n_cropped + 1
  }
  cat("Cropped", n_cropped, "of", length(ppt_files), "ppt file(s) to bbox_extent.\n\n")
}

crop_ppt_archive_inplace()

# ============================================================
# 4. Check for the Kling fog time series raster
#    Source data ships as Kling/timeseries.tif.zip and must be decompressed
#    once, outside this script, e.g.:
#      cd /private/groups/meyerlab/esteiger/env/fog/Kling && \
#        unzip -o timeseries.tif.zip timeseries.tif
# ============================================================

if (!file.exists(kling_ts_file)) {
  stop(
    "Kling fog time series raster not found at ", kling_ts_file, ".\n",
    "Decompress it first, e.g.:\n",
    "  cd ", dirname(kling_ts_file), " && unzip -o ", basename(kling_ts_zip),
    " ", basename(kling_ts_file)
  )
}

r_kling_ts  <- rast(kling_ts_file)
kling_years <- sort(as.integer(names(r_kling_ts)))
cat("Kling fog time series:", length(kling_years), "years (",
    min(kling_years), "-", max(kling_years), ")\n\n")

# ============================================================
# 5. Helper functions
# ============================================================

# --- fog: mean FLCC across Jun-Sep for a single year ---
extract_fog_seasonal <- function(year, lats, lons) {
  paths <- file.path(fog_dir, sprintf("%d-%02d.tif", year, fog_months))
  names(paths) <- sprintf("%d-%02d", year, fog_months)

  missing_months <- names(paths)[!file.exists(paths)]
  if (length(missing_months) > 0) {
    stop("Missing fog raster(s) for year ", year, ": ",
         paste(missing_months, collapse = ", "))
  }

  pts <- vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326")

  monthly_vals <- matrix(
    sapply(unname(paths), function(p) {
      r        <- rast(p)
      pts_proj <- project(pts, crs(r))
      extract(r, pts_proj)[, 2]
    }),
    nrow = length(lats)
  )

  rowMeans(monthly_vals, na.rm = FALSE)
}

# --- fog: matrix of annual values across multiple years (cols = years) ---
extract_fog_multi_year <- function(years, lats, lons) {
  n    <- length(lats)
  vals <- matrix(NA_real_, nrow = n, ncol = length(years), dimnames = list(NULL, years))
  for (i in seq_along(years)) {
    vals[, i] <- tryCatch(extract_fog_seasonal(years[i], lats, lons), error = function(e) {
      warning("Could not extract fog data for year ", years[i], ": ", e$message)
      rep(NA_real_, n)
    })
  }
  vals
}

# --- PRISM: BIO12/15/18 for a single water year ---
extract_prism_year <- function(sample_year, lats, lons) {
  wy_yyyymm <- c(
    sprintf("%04d%02d", sample_year - 1, 10:12),
    sprintf("%04d%02d", sample_year,     1:9)
  )
  bio18_yyyymm <- sprintf("%04d%02d", sample_year, 6:8)  # Jun–Aug of sampling year

  bil_paths <- sapply(wy_yyyymm, find_prism_bil)
  missing   <- wy_yyyymm[sapply(bil_paths, is.null)]
  if (length(missing) > 0) {
    stop("Missing PRISM raster(s) for water year ", sample_year, ": ",
         paste(missing, collapse = ", "))
  }

  stack        <- rast(unlist(bil_paths))
  names(stack) <- wy_yyyymm

  pts <- vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326")
  pts <- project(pts, crs(stack))

  vals <- as.matrix(extract(stack, pts)[, -1])

  bio12 <- rowSums(vals)
  bio15 <- apply(vals, 1, function(x) {
    if (any(is.na(x))) return(NA_real_)
    (sd(x) / mean(x)) * 100
  })
  bio18 <- rowSums(vals[, bio18_yyyymm, drop = FALSE])

  list(BIO12 = bio12, BIO15 = bio15, BIO18 = bio18)
}

# --- PRISM: BIO12/15/18 across multiple water years (cols = years) ---
extract_prism_multi_year <- function(years, lats, lons) {
  n <- length(lats)
  bio12 <- bio15 <- bio18 <- matrix(NA_real_, nrow = n, ncol = length(years),
                                     dimnames = list(NULL, years))
  for (i in seq_along(years)) {
    pv <- tryCatch(extract_prism_year(years[i], lats, lons), error = function(e) {
      warning("Could not extract PRISM data for water year ", years[i], ": ", e$message)
      NULL
    })
    if (!is.null(pv)) {
      bio12[, i] <- pv$BIO12
      bio15[, i] <- pv$BIO15
      bio18[, i] <- pv$BIO18
    }
  }
  list(BIO12 = bio12, BIO15 = bio15, BIO18 = bio18)
}

# --- Kling: annual fog values from the real time series (cols = years) ---
extract_kling_years <- function(years, lats, lons) {
  n       <- length(lats)
  yrs_chr <- as.character(years)
  vals    <- matrix(NA_real_, nrow = n, ncol = length(years), dimnames = list(NULL, yrs_chr))

  present <- yrs_chr %in% names(r_kling_ts)
  if (any(present)) {
    pts <- vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326")
    pts <- project(pts, crs(r_kling_ts))
    ex  <- as.matrix(extract(r_kling_ts[[yrs_chr[present]]], pts)[, -1, drop = FALSE])
    vals[, present] <- ex
  }
  vals
}

# ============================================================
# 5b. Environmental variable selection: bioclim analogs, collinearity, and
#     PCA — all at the _25yr timescale only (not _germyr or _<N>yr).
#
#     Adds analogs of the remaining WorldClim/Bioclim variables (beyond
#     BIO12/15/18, already computed above) from the same PRISM data:
#       BIO1  Annual mean temperature          BIO9  Mean temp, driest qtr
#       BIO2  Mean diurnal range               BIO10 Mean temp, warmest qtr
#       BIO3  Isothermality (BIO2/BIO7*100)    BIO11 Mean temp, coldest qtr
#       BIO4  Temperature seasonality (SD*100) BIO13 Precip, wettest month
#       BIO5  Max temp, warmest month          BIO14 Precip, driest month
#       BIO6  Min temp, coldest month          BIO16 Precip, wettest qtr
#       BIO7  Temperature annual range         BIO17 Precip, driest qtr
#       BIO8  Mean temp, wettest qtr           BIO19 Precip, coldest qtr
#     "Quarter" = a rolling 3-consecutive-month window (12 such windows per
#     year, circular), the standard BIOCLIM definition — not the fixed
#     Jun-Aug window used for the existing BIO18.
#     Requires PRISM tmean/tmax/tmin (not otherwise used in this script);
#     downloaded pre-cropped to bbox_extent, for the _25yr reference window
#     only (water years 2000-2024).
#
#     NOTE: this section writes diagnostic plots (heatmap/scree/biplot) but
#     no txt/csv — see "remove_collinear_25yr" below, to be populated after
#     reviewing the flagged pairs together.
# ============================================================

cat("\n=== Section 5b: bioclim analogs + collinearity/PCA at _25yr ===\n\n")

# --- 5b.1 Download + crop PRISM tmean/tmax/tmin for the _25yr window ---

temp_vars <- c("tmean", "tmax", "tmin")

wy_months_25yr <- do.call(rbind, lapply(ref_bio_years, function(yr) {
  rbind(
    data.frame(year = yr - 1, month = 10:12),
    data.frame(year = yr,     month = 1:9)
  )
}))
wy_months_25yr$yyyymm <- sprintf("%04d%02d", wy_months_25yr$year, wy_months_25yr$month)
unique_yyyymm_25yr    <- sort(unique(wy_months_25yr$yyyymm))

find_prism_bil_var <- function(var, yyyymm) {
  dir_path <- file.path(prism_dir, yyyymm)
  files <- list.files(dir_path, pattern = paste0("^prism_", var, "_us_30s_.*\\.(tif|bil)$"),
                       full.names = TRUE)
  if (length(files) == 0) return(NULL)
  files[1]
}

# PRISM's web service enforces a hard per-file limit ("more than twice in
# one day (Pacific local time)... repeated offenses may result in your IP
# address being blocked") and can also return transient 500s under load.
# These need different handling: a rate-limit response is NOT retriable
# today (retrying just adds more offenses against the same file and risks
# an IP block for the whole shared server) and aborts the whole download
# pass immediately; a genuine transient error gets a couple of gentle
# retries. Requests are paced with Sys.sleep() between files either way.
prism_rate_limited <- FALSE

download_and_crop_prism_var <- function(var, yyyymm, max_tries = 2) {
  url <- paste0("https://services.nacse.org/prism/data/get/us/800m/", var, "/", yyyymm)
  for (attempt in seq_len(max_tries)) {
    zip_path <- tempfile(fileext = ".zip")
    result <- tryCatch({
      download.file(url, zip_path, mode = "wb", quiet = TRUE)
      first_bytes <- readBin(zip_path, "raw", n = 4)
      is_zip <- length(first_bytes) == 4 && identical(first_bytes[1:2], as.raw(c(0x50, 0x4b)))
      if (!is_zip) {
        body <- tryCatch(readLines(zip_path, n = 5, warn = FALSE), error = function(e) "")
        if (any(grepl("more than twice|IP address being blocked", body, ignore.case = TRUE))) {
          return(list(status = "rate_limited", body = paste(body, collapse = " ")))
        }
        stop("response was not a zip file: ", paste(body, collapse = " "))
      }
      tmp_dir <- tempfile()
      dir.create(tmp_dir)
      unzip(zip_path, exdir = tmp_dir)
      unlink(zip_path)
      full_tif <- list.files(tmp_dir, pattern = "_30s_.*\\.(tif|bil)$", full.names = TRUE)[1]
      if (is.na(full_tif)) stop("no raster found in downloaded zip")
      r <- crop(rast(full_tif), bbox_extent)
      out_dir <- file.path(prism_dir, yyyymm)
      dir.create(out_dir, showWarnings = FALSE)
      out_path <- file.path(out_dir, sprintf("prism_%s_us_30s_%s.tif", var, yyyymm))
      writeRaster(r, out_path, overwrite = TRUE)
      unlink(tmp_dir, recursive = TRUE)
      list(status = "ok")
    }, error = function(e) {
      unlink(zip_path)
      list(status = "error", body = e$message)
    })
    if (identical(result$status, "ok")) return(invisible("ok"))
    if (identical(result$status, "rate_limited")) {
      warning("PRISM rate limit hit for ", var, " ", yyyymm, ": ", result$body)
      return(invisible("rate_limited"))
    }
    if (attempt < max_tries) {
      cat("    retry for", var, yyyymm, "after transient error:", result$body, "\n")
      Sys.sleep(3)
    }
  }
  warning("Failed to download/crop PRISM ", var, " for ", yyyymm,
          " after ", max_tries, " attempts.")
  invisible("error")
}

for (var in temp_vars) {
  if (prism_rate_limited) break
  cat("Checking PRISM", var, "files (25yr reference window)...\n")
  missing_var <- unique_yyyymm_25yr[sapply(unique_yyyymm_25yr,
                                            function(m) is.null(find_prism_bil_var(var, m)))]
  if (length(missing_var) == 0) {
    cat("  All required", var, "files already present.\n\n")
  } else {
    cat("  Downloading", length(missing_var), "missing", var, "file(s)...\n")
    for (yyyymm in missing_var) {
      cat("    Downloading:", var, yyyymm, "\n")
      status <- download_and_crop_prism_var(var, yyyymm)
      if (identical(status, "rate_limited")) {
        cat("  PRISM rate limit reached — stopping further downloads for today.\n",
            "  Re-run the script another day; already-downloaded files are cached",
            " and won't be re-requested.\n")
        prism_rate_limited <- TRUE
        break
      }
      Sys.sleep(1)
    }
    cat("\n")
  }
}

# --- 5b.2 Generic monthly-matrix extraction + rolling-quarter bioclim analogs ---

extract_prism_monthly_matrix <- function(var, sample_year, lats, lons) {
  wy_yyyymm <- c(sprintf("%04d%02d", sample_year - 1, 10:12),
                 sprintf("%04d%02d", sample_year,     1:9))
  bil_paths <- sapply(wy_yyyymm, function(m) find_prism_bil_var(var, m))
  missing   <- wy_yyyymm[sapply(bil_paths, is.null)]
  if (length(missing) > 0) {
    stop("Missing PRISM ", var, " raster(s) for water year ", sample_year, ": ",
         paste(missing, collapse = ", "))
  }
  stack        <- rast(unlist(bil_paths))
  names(stack) <- wy_yyyymm
  pts <- vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326")
  pts <- project(pts, crs(stack))
  as.matrix(extract(stack, pts)[, -1])
}

# Rolling 3-month window aggregation, circular over the 12 water-year months
# (valid regardless of which month the vector starts on).
roll3 <- function(mat, fun) {
  n <- ncol(mat)
  sapply(seq_len(n), function(i) {
    idx <- ((i - 1):(i + 1)) %% n + 1
    apply(mat[, idx, drop = FALSE], 1, fun)
  })
}

# Bioclim analogs for a single water year (BIO12/15/18 excluded — already
# computed by extract_prism_year() above).
compute_bioclim_analogs <- function(sample_year, lats, lons) {
  ppt   <- extract_prism_monthly_matrix("ppt",   sample_year, lats, lons)
  tmean <- extract_prism_monthly_matrix("tmean", sample_year, lats, lons)
  tmax  <- extract_prism_monthly_matrix("tmax",  sample_year, lats, lons)
  tmin  <- extract_prism_monthly_matrix("tmin",  sample_year, lats, lons)

  bio1 <- rowMeans(tmean)
  bio2 <- rowMeans(tmax - tmin)
  bio5 <- apply(tmax, 1, max)
  bio6 <- apply(tmin, 1, min)
  bio7 <- bio5 - bio6
  bio3 <- ifelse(bio7 == 0, NA_real_, (bio2 / bio7) * 100)
  bio4 <- apply(tmean, 1, sd) * 100
  bio13 <- apply(ppt, 1, max)
  bio14 <- apply(ppt, 1, min)

  qtr_temp   <- roll3(tmean, mean)  # n x 12: mean temp of each rolling quarter
  qtr_precip <- roll3(ppt,   sum)   # n x 12: total precip of each rolling quarter

  bio10 <- apply(qtr_temp, 1, max)
  bio11 <- apply(qtr_temp, 1, min)
  bio16 <- apply(qtr_precip, 1, max)
  bio17 <- apply(qtr_precip, 1, min)
  n_row <- nrow(qtr_temp)
  bio8  <- qtr_temp[cbind(seq_len(n_row), max.col(qtr_precip,  ties.method = "first"))]
  bio9  <- qtr_temp[cbind(seq_len(n_row), max.col(-qtr_precip, ties.method = "first"))]
  bio19 <- qtr_precip[cbind(seq_len(n_row), max.col(-qtr_temp, ties.method = "first"))]

  list(BIO1 = bio1, BIO2 = bio2, BIO3 = bio3, BIO4 = bio4, BIO5 = bio5, BIO6 = bio6,
       BIO7 = bio7, BIO8 = bio8, BIO9 = bio9, BIO10 = bio10, BIO11 = bio11,
       BIO13 = bio13, BIO14 = bio14, BIO16 = bio16, BIO17 = bio17, BIO19 = bio19)
}

compute_bioclim_25yr <- function(years, lats, lons) {
  n <- length(lats)
  bio_names <- c("BIO1", "BIO2", "BIO3", "BIO4", "BIO5", "BIO6", "BIO7", "BIO8", "BIO9",
                 "BIO10", "BIO11", "BIO13", "BIO14", "BIO16", "BIO17", "BIO19")
  acc <- setNames(vector("list", length(bio_names)), bio_names)
  for (b in bio_names) acc[[b]] <- matrix(NA_real_, nrow = n, ncol = length(years))
  for (i in seq_along(years)) {
    cat("  Bioclim analogs, water year", years[i], "...\n")
    bv <- tryCatch(compute_bioclim_analogs(years[i], lats, lons), error = function(e) {
      warning("Could not compute bioclim analogs for water year ", years[i], ": ", e$message)
      NULL
    })
    if (!is.null(bv)) for (b in bio_names) acc[[b]][, i] <- bv[[b]]
  }
  out <- sapply(bio_names, function(b) rowMeans(acc[[b]], na.rm = TRUE))
  colnames(out) <- paste0(bio_names, "_25yr")
  out
}

# --- 5b.3 Locality-equalized coordinate set, then assemble the _25yr
#     analysis table (one row per unique coordinate)
#     Samples/coordinates/localities come from gvtp_phenotype_13feb2026.tsv
#     (the phenotype dataset, via its "Locality" column), not
#     all_coord_gvtp.csv, per user request. Each locality is capped at
#     MAX_PER_LOCALITY samples so no single locality (e.g. Strauss n=50,
#     VAFB_Lion_head n=48) dominates the collinearity/PCA; localities at or
#     below the cap (e.g. Dangermond n=5) are kept in full rather than
#     downsampled further, so they aren't drastically under-represented.
# fog_25yr / BIO12_25yr / BIO15_25yr / BIO18_25yr / CHILI are recomputed here
# (same helper functions/inputs as section 7 below) so this selection
# analysis is self-contained near the top of the script, as requested.
# Kling has no "_25yr" column (see header) and is not included here.

MAX_PER_LOCALITY <- 15

pheno_sel <- read.delim(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)
pheno_sel <- pheno_sel[, c("Locality", "Lat _dd.dddd", "Lon_dd.dddd")]
colnames(pheno_sel) <- c("Locality", "Lat", "Long")
pheno_sel <- pheno_sel[!is.na(pheno_sel$Locality) & pheno_sel$Locality != "" &
                       !is.na(pheno_sel$Lat) & !is.na(pheno_sel$Long), ]

cat("\nSelection analysis: locality sample sizes before capping at",
    MAX_PER_LOCALITY, ":\n")
print(table(pheno_sel$Locality))

set.seed(42)
capped_sel <- do.call(rbind, lapply(split(pheno_sel, pheno_sel$Locality), function(d) {
  if (nrow(d) > MAX_PER_LOCALITY) d[sample(nrow(d), MAX_PER_LOCALITY), ] else d
}))
rownames(capped_sel) <- NULL

cat("\nLocality sample sizes after capping:\n")
print(table(capped_sel$Locality))

uniq_xy_sel <- unique(capped_sel[, c("Lat", "Long")])
rownames(uniq_xy_sel) <- NULL
cat("\nSelection analysis: ", nrow(uniq_xy_sel),
    "unique coordinates (from", nrow(capped_sel), "locality-capped samples)\n")

cat("  fog_25yr...\n")
fog_mat_sel <- extract_fog_multi_year(ref_fog_years, uniq_xy_sel$Lat, uniq_xy_sel$Long)
fog_25yr_sel <- rowMeans(fog_mat_sel, na.rm = TRUE)

cat("  BIO12/15/18_25yr...\n")
bio_mat_sel <- extract_prism_multi_year(ref_bio_years, uniq_xy_sel$Lat, uniq_xy_sel$Long)
ref_cols_sel <- as.character(ref_bio_years)
BIO12_25yr_sel <- rowMeans(bio_mat_sel$BIO12[, ref_cols_sel], na.rm = TRUE)
BIO15_25yr_sel <- rowMeans(bio_mat_sel$BIO15[, ref_cols_sel], na.rm = TRUE)
BIO18_25yr_sel <- rowMeans(bio_mat_sel$BIO18[, ref_cols_sel], na.rm = TRUE)

cat("  CHILI...\n")
chili_rast_sel <- rast(chili_file)
pts_sel <- vect(data.frame(x = uniq_xy_sel$Long, y = uniq_xy_sel$Lat),
                geom = c("x", "y"), crs = "EPSG:4326")
CHILI_sel <- extract(chili_rast_sel, project(pts_sel, crs(chili_rast_sel)))[, 2]

cat("  BIO1-11/13/14/16/17/19_25yr (new analogs)...\n")
bioclim_25yr_sel <- compute_bioclim_25yr(ref_bio_years, uniq_xy_sel$Lat, uniq_xy_sel$Long)

sel_vars <- data.frame(
  Lat  = uniq_xy_sel$Lat,
  Long = uniq_xy_sel$Long,
  fog_25yr   = fog_25yr_sel,
  CHILI      = CHILI_sel,
  BIO12_25yr = BIO12_25yr_sel,
  BIO15_25yr = BIO15_25yr_sel,
  BIO18_25yr = BIO18_25yr_sel,
  bioclim_25yr_sel,
  check.names = FALSE
)

# --- 5b.4 Collinearity — Spearman correlations among all _25yr variables ---

COLLINEAR_THRESH_25yr <- 0.7
sel_num_vars <- setdiff(names(sel_vars), c("Lat", "Long"))

cor_mat_sel25 <- cor(sel_vars[, sel_num_vars], method = "spearman",
                     use = "pairwise.complete.obs")

# Variables with no variance (or insufficient paired observations, e.g. a
# failed download leaving a column all-NA) produce an all-NA row/column in
# cor_mat_sel25, which crashes hclust below — drop them here rather than
# silently including undefined correlations.
all_na_row    <- apply(cor_mat_sel25, 1, function(x) sum(!is.na(x)) <= 1)  # only the diagonal is non-NA
na_vars_sel25 <- rownames(cor_mat_sel25)[all_na_row]
if (length(na_vars_sel25) > 0) {
  cat("\nDropping variable(s) with no valid pairwise correlations",
      "(all-NA/constant — check upstream extraction):\n")
  cat(" ", paste(na_vars_sel25, collapse = ", "), "\n")
  keep_sel25    <- setdiff(rownames(cor_mat_sel25), na_vars_sel25)
  cor_mat_sel25 <- cor_mat_sel25[keep_sel25, keep_sel25]
  sel_num_vars  <- keep_sel25
}

pairs_idx_sel25 <- which(upper.tri(cor_mat_sel25), arr.ind = TRUE)
cor_pairs_sel25 <- data.frame(
  var1 = rownames(cor_mat_sel25)[pairs_idx_sel25[, 1]],
  var2 = colnames(cor_mat_sel25)[pairs_idx_sel25[, 2]],
  rho  = round(cor_mat_sel25[pairs_idx_sel25], 3)
)
cor_pairs_sel25 <- cor_pairs_sel25[order(-abs(cor_pairs_sel25$rho)), ]
flagged_sel25   <- cor_pairs_sel25[!is.na(cor_pairs_sel25$rho) &
                                    abs(cor_pairs_sel25$rho) >= COLLINEAR_THRESH_25yr, ]

cat(sprintf("\n=== _25yr collinearity (Spearman |rho| >= %.1f) ===\n", COLLINEAR_THRESH_25yr))
cat(sprintf("Flagged pairs: %d\n\n", nrow(flagged_sel25)))
print(flagged_sel25, row.names = FALSE)

hc_ord_sel25     <- hclust(as.dist(1 - cor_mat_sel25))$order
var_levels_sel25 <- rownames(cor_mat_sel25)[hc_ord_sel25]

cor_df_sel25 <- data.frame(
  var1 = rep(rownames(cor_mat_sel25), ncol(cor_mat_sel25)),
  var2 = rep(colnames(cor_mat_sel25), each = nrow(cor_mat_sel25)),
  rho  = as.vector(cor_mat_sel25)
)
cor_df_sel25$var1 <- factor(cor_df_sel25$var1, levels = var_levels_sel25)
cor_df_sel25$var2 <- factor(cor_df_sel25$var2, levels = var_levels_sel25)

jpeg(file.path(plots_dir, "env_selection_correlation_heatmap_25yr.jpg"),
     width = 1600, height = 1400)
print(
  ggplot(cor_df_sel25, aes(x = var2, y = var1, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title   = element_blank()) +
    ggtitle("_25yr env variable Spearman correlations (fog/CHILI/BIO12-19 analogs)")
)
dev.off()
cat(sprintf("Heatmap written to %s\n",
            file.path(plots_dir, "env_selection_correlation_heatmap_25yr.jpg")))

# --- 5b.5 PCA on the retained _25yr variables ---
#     remove_collinear_25yr decisions so far (from a precip-only preview
#     pass, locality-capped subsample, n=63 unique coordinates — see
#     scratch script precip_only_25yr_collinearity.R; BIO1-11 pending the
#     PRISM temp downloads and not yet reviewed):
#
#     Wet-season cluster — BIO12/BIO13/BIO16/BIO19 all rho ~= 0.99 with each
#     other. In this Mediterranean climate nearly all precip falls in
#     winter, so annual total (BIO12), wettest month (BIO13), wettest
#     quarter (BIO16), and the fixed winter quarter (BIO19) are all
#     measuring essentially the same signal.
#       Kept:    BIO12_25yr (annual precip — most standard/interpretable)
#       Removed: BIO13_25yr (rho = 0.991 with BIO12_25yr)
#                BIO16_25yr (rho = 0.994 with BIO12_25yr)
#                BIO19_25yr (rho = 0.988 with BIO12_25yr)
#
#     Dry-season cluster — BIO14/BIO17/BIO18 rho = 0.84-0.94 with each
#     other (dry season ~= summer here, so driest month/quarter and the
#     fixed summer quarter all track together).
#       Kept:    BIO18_25yr (precip of warmest/summer quarter — user pick)
#       Removed: BIO14_25yr (rho = 0.840 with BIO18_25yr)
#                BIO17_25yr (rho = 0.861 with BIO18_25yr)
#
#     BIO15_25yr (seasonality CV) stayed under the 0.7 threshold (max
#     rho = 0.68, with the wet cluster) and is kept. fog_25yr and CHILI
#     were not strongly collinear with anything in the precip-only pass.
#
#     Re-review and extend this list once the full run (including BIO1-11)
#     completes.

remove_collinear_25yr <- c(
  "BIO13_25yr",  # wet-season cluster — remove instead of BIO12_25yr (rho = 0.991)
  "BIO16_25yr",  # wet-season cluster — remove instead of BIO12_25yr (rho = 0.994)
  "BIO19_25yr",  # wet-season cluster — remove instead of BIO12_25yr (rho = 0.988)
  "BIO14_25yr",  # dry-season cluster — remove instead of BIO18_25yr (rho = 0.840)
  "BIO17_25yr"   # dry-season cluster — remove instead of BIO18_25yr (rho = 0.861)
)

pca_vars_sel25 <- setdiff(sel_num_vars, remove_collinear_25yr)
pca_data_sel25 <- sel_vars[, pca_vars_sel25]
pca_data_sel25 <- pca_data_sel25[complete.cases(pca_data_sel25), ]
cat(sprintf("\n_25yr PCA: %d variables, %d complete-case coordinates\n",
            ncol(pca_data_sel25), nrow(pca_data_sel25)))

pca_res_sel25 <- prcomp(pca_data_sel25, center = TRUE, scale. = TRUE)

var_exp_sel25 <- summary(pca_res_sel25)$importance
cat("\n_25yr variance explained:\n")
print(round(var_exp_sel25[, 1:min(10, ncol(var_exp_sel25))], 3))

n_pc_sel25 <- min(10, ncol(pca_res_sel25$rotation))
cat("\n_25yr loadings (PC1-PC", n_pc_sel25, "):\n", sep = "")
print(round(pca_res_sel25$rotation[, 1:n_pc_sel25], 3))

pct_var_sel25  <- 100 * pca_res_sel25$sdev^2 / sum(pca_res_sel25$sdev^2)
scree_df_sel25 <- data.frame(PC = seq_along(pct_var_sel25), pct_var = pct_var_sel25,
                             cumvar = cumsum(pct_var_sel25))

jpeg(file.path(plots_dir, "env_selection_pca_scree_25yr.jpg"), width = 900, height = 700)
print(
  ggplot(scree_df_sel25, aes(x = PC, y = pct_var)) +
    geom_col(fill = "steelblue") +
    geom_line(aes(y = cumvar)) +
    geom_point(aes(y = cumvar)) +
    geom_hline(yintercept = 80, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = seq_len(nrow(scree_df_sel25))) +
    labs(x = "PC", y = "Variance explained (%)") +
    theme_bw(base_size = 14)
)
dev.off()
cat(sprintf("Scree plot written to %s\n",
            file.path(plots_dir, "env_selection_pca_scree_25yr.jpg")))

scores_sel25 <- as.data.frame(pca_res_sel25$x)
arrows_df_sel25 <- as.data.frame(pca_res_sel25$rotation[, 1:2])
arrows_df_sel25$var <- rownames(arrows_df_sel25)
scale_arrow_sel25 <- 4

jpeg(file.path(plots_dir, "env_selection_pca_biplot_25yr.jpg"), width = 1100, height = 1000)
print(
  ggplot(scores_sel25, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.6, size = 2, colour = "steelblue") +
    geom_segment(data = arrows_df_sel25,
                 aes(x = 0, y = 0,
                     xend = PC1 * scale_arrow_sel25, yend = PC2 * scale_arrow_sel25),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey30") +
    geom_text(data = arrows_df_sel25,
              aes(x = PC1 * scale_arrow_sel25 * 1.1, y = PC2 * scale_arrow_sel25 * 1.1,
                  label = var),
              inherit.aes = FALSE, size = 4, colour = "grey20") +
    labs(x = sprintf("PC1 (%.1f%%)", pct_var_sel25[1]),
         y = sprintf("PC2 (%.1f%%)", pct_var_sel25[2])) +
    theme_bw(base_size = 14)
)
dev.off()
cat(sprintf("Biplot written to %s\n\n",
            file.path(plots_dir, "env_selection_pca_biplot_25yr.jpg")))

# ============================================================
# 6. Per-sample "_germyr" values: fog, CHILI, BIO12/15/18, Kling for the
#    specific year each plant was collected
# ============================================================

cat("Extracting germination-year (per-sample) values...\n")

results <- coords
results$fog_germyr   <- NA_real_
results$BIO12_germyr <- NA_real_
results$BIO15_germyr <- NA_real_
results$BIO18_germyr <- NA_real_
results$Kling_germyr <- NA_real_

for (yr in sample_years) {
  idx  <- which(results$Year == yr)
  lats <- results$Lat[idx]
  lons <- results$Long[idx]
  cat("  Year", yr, "—", length(idx), "sample(s)\n")

  if (yr >= fog_year_min && yr <= fog_year_max) {
    results$fog_germyr[idx] <- extract_fog_seasonal(yr, lats, lons)
  }

  pv <- tryCatch(extract_prism_year(yr, lats, lons), error = function(e) {
    warning("Could not extract PRISM data for year ", yr, ": ", e$message)
    NULL
  })
  if (!is.null(pv)) {
    results$BIO12_germyr[idx] <- pv$BIO12
    results$BIO15_germyr[idx] <- pv$BIO15
    results$BIO18_germyr[idx] <- pv$BIO18
  }

  results$Kling_germyr[idx] <- extract_kling_years(yr, lats, lons)[, 1]
}

# --- CHILI: static layer, single column, no time variants ---

cat("\nExtracting CHILI values...\n")

chili_rast <- rast(chili_file)
pts_all <- vect(
  data.frame(x = results$Long, y = results$Lat),
  geom = c("x", "y"), crs = "EPSG:4326"
)
pts_chili     <- project(pts_all, crs(chili_rast))
results$CHILI <- extract(chili_rast, pts_chili)[, 2]

cat("  CHILI extracted for", sum(!is.na(results$CHILI)), "/", nrow(results), "samples\n")

# ============================================================
# 7. Fixed-window multi-year averages ("_25yr", "_<N>yr")
#    Computed once per unique coordinate, then joined back onto every
#    sample sharing that coordinate. The "_<N>yr" long window is sized
#    and dated to match Kling fog's real observed range (kling_years,
#    loaded in section 4) exactly, so BIO12/15/18_<N>yr and Kling_<N>yr
#    are averages over the identical set of years, not just the same
#    N — the previous version averaged BIO12/15/18 over the full
#    1927-2024 PRISM archive (N=98) while Kling only had 1950-2014
#    (N=65), so the two "long" windows didn't actually overlap in years.
# ============================================================

cat("\nExtracting fixed-window multi-year averages (per unique coordinate)...\n")

uniq_xy <- unique(results[, c("Lat", "Long")])
rownames(uniq_xy) <- NULL
cat("  Unique coordinates:", nrow(uniq_xy), "\n")

# --- fog_25yr: mean of annual Jun-Sep fog, 2000-2025 ---
cat("  fog_25yr (", min(ref_fog_years), "-", max(ref_fog_years), ")...\n")
fog_mat <- extract_fog_multi_year(ref_fog_years, uniq_xy$Lat, uniq_xy$Long)
uniq_xy$fog_25yr <- rowMeans(fog_mat, na.rm = TRUE)

# --- BIO12/15/18: one pass over the years needed for both the fixed
#     "_25yr" window (2000-2024) and the Kling-matched long window
#     (kling_years), from which both are derived. PRISM files for the
#     full 1927-2024 archive stay on disk (still used by section 10) —
#     this just limits which cached years get read/averaged here. ---
bio_years_needed <- sort(union(ref_bio_years, kling_years))

cat("  BIO12/15/18 (", min(bio_years_needed), "-", max(bio_years_needed), ")...\n")
bio_mat   <- extract_prism_multi_year(bio_years_needed, uniq_xy$Lat, uniq_xy$Long)
ref_cols  <- as.character(ref_bio_years)
long_cols <- as.character(kling_years)

uniq_xy[["BIO12_25yr"]] <- rowMeans(bio_mat$BIO12[, ref_cols], na.rm = TRUE)
uniq_xy[["BIO15_25yr"]] <- rowMeans(bio_mat$BIO15[, ref_cols], na.rm = TRUE)
uniq_xy[["BIO18_25yr"]] <- rowMeans(bio_mat$BIO18[, ref_cols], na.rm = TRUE)

# Long-window suffix/years now match Kling exactly (e.g. "65yr" over
# 1950-2014), replacing the old mismatched "98yr" over 1927-2024.
bio_all_suffix <- paste0(length(kling_years), "yr")
uniq_xy[[paste0("BIO12_", bio_all_suffix)]] <- rowMeans(bio_mat$BIO12[, long_cols], na.rm = TRUE)
uniq_xy[[paste0("BIO15_", bio_all_suffix)]] <- rowMeans(bio_mat$BIO15[, long_cols], na.rm = TRUE)
uniq_xy[[paste0("BIO18_", bio_all_suffix)]] <- rowMeans(bio_mat$BIO18[, long_cols], na.rm = TRUE)

# --- Kling: mean over the full observed time series — same years as
#     BIO12/15/18_<N>yr above, so the suffix matches too ---
kling_all_suffix <- bio_all_suffix

cat("  Kling (", min(kling_years), "-", max(kling_years), ")...\n")
kling_mat <- extract_kling_years(kling_years, uniq_xy$Lat, uniq_xy$Long)
uniq_xy[[paste0("Kling_", kling_all_suffix)]] <- rowMeans(kling_mat, na.rm = TRUE)

# --- join back onto every sample sharing that coordinate ---
results <- merge(results, uniq_xy, by = c("Lat", "Long"), all.x = TRUE, sort = FALSE)
results <- results[match(coords$SampleID, results$SampleID), ]
rownames(results) <- NULL

# ============================================================
# 8. Diagnostic plots — coordinates overlaid on each env layer
# ============================================================

cat("\nGenerating diagnostic plots...\n")

pts_plot <- vect(
  data.frame(x = results$Long, y = results$Lat),
  geom = c("x", "y"), crs = "EPSG:4326"
)

# --- Fog: one representative layer (2020-06) ---
fog_check_path <- file.path(fog_dir, "2020-06.tif")
if (file.exists(fog_check_path)) {
  png(file.path(plots_dir, "check_fog_coords.png"), width = 1200, height = 900, res = 150)
  fog_check <- rast(fog_check_path)
  plot(fog_check, main = "Fog check: 2020-06 with sample coordinates",
       col = hcl.colors(100, "Blues 3", rev = TRUE))
  points(project(pts_plot, crs(fog_check)), pch = 21, bg = "red", col = "black", cex = 0.8)
  dev.off()
  cat("  Saved: check_fog_coords.png\n")
} else {
  warning("Could not find 2020-06.tif for fog diagnostic plot.")
}

# --- CHILI ---
png(file.path(plots_dir, "check_chili_coords.png"), width = 1200, height = 900, res = 150)
plot(chili_rast, main = "CHILI check: sample coordinates",
     col = hcl.colors(100, "RdYlGn"))
points(project(pts_plot, crs(chili_rast)), pch = 21, bg = "red", col = "black", cex = 0.8)
dev.off()
cat("  Saved: check_chili_coords.png\n")

# --- PRISM: one representative monthly raster ---
prism_check_yyyymm <- "202006"
prism_check_bil    <- find_prism_bil(prism_check_yyyymm)
if (!is.null(prism_check_bil)) {
  png(file.path(plots_dir, "check_prism_coords.png"), width = 1200, height = 900, res = 150)
  prism_check <- rast(prism_check_bil)
  plot(prism_check, main = "PRISM check: Jun 2020 precipitation with sample coordinates",
       col = hcl.colors(100, "Blues", rev = TRUE))
  points(project(pts_plot, crs(prism_check)), pch = 21, bg = "red", col = "black", cex = 0.8)
  dev.off()
  cat("  Saved: check_prism_coords.png\n")
} else {
  warning("Could not find PRISM 202006 raster for diagnostic plot.")
}

# --- Kling: one representative annual layer ---
kling_check_year <- as.character(kling_years[which.min(abs(kling_years - 2000))])
png(file.path(plots_dir, "check_kling_coords.png"), width = 1200, height = 900, res = 150)
plot(r_kling_ts[[kling_check_year]],
     main = paste("Kling fog check:", kling_check_year, "with sample coordinates"),
     col = hcl.colors(100, "Blues 3", rev = TRUE))
points(project(pts_plot, crs(r_kling_ts)), pch = 21, bg = "red", col = "black", cex = 0.8)
dev.off()
cat("  Saved: check_kling_coords.png\n")

# ============================================================
# 9. Save output
# ============================================================

output <- results[, c(
  "SampleID", "Year", "Lat", "Long",
  "fog_germyr", "fog_25yr",
  "CHILI",
  "BIO12_germyr", "BIO12_25yr", paste0("BIO12_", bio_all_suffix),
  "BIO15_germyr", "BIO15_25yr", paste0("BIO15_", bio_all_suffix),
  "BIO18_germyr", "BIO18_25yr", paste0("BIO18_", bio_all_suffix),
  "Kling_germyr", paste0("Kling_", kling_all_suffix)
)]
write.csv(output, output_file, row.names = FALSE)

cat("\nResults saved to:", output_file, "\n")
for (col in setdiff(names(output), c("SampleID", "Year", "Lat", "Long"))) {
  cat(sprintf("  %-16s %d / %d\n", col, sum(!is.na(output[[col]])), nrow(output)))
}
if (length(years_out_range) > 0) {
  cat("Samples with NA fog_germyr (year outside fog dataset range):",
      sum(output$Year %in% years_out_range), "\n")
}

# ============================================================
# 10. Extract water-year PRISM + fog variables across all years
#     for every unique tarplant coordinate
#     Input:  gvtp_phenotype_13feb2026.tsv  (Lat _dd.dddd, Lon_dd.dddd, Locality)
#     Output: env_wi_range.csv
#             Columns: Locality, Lat, Long,
#                      BIO12_1927 … BIO12_<last_yr>,
#                      BIO15_1927 … BIO15_<last_yr>,
#                      BIO18_1927 … BIO18_<last_yr>,
#                      fog_2000 … fog_2025,
#                      kling_fog_1950 … kling_fog_2014
# ============================================================

range_output <- "/private/groups/meyerlab/esteiger/pheno/env_wi_range.csv"

# --- helpers (redefined here so this section runs independently) ---

find_prism_bil_r9 <- function(yyyymm) {
  dir_path <- file.path(prism_dir, yyyymm)
  files    <- list.files(dir_path, pattern = "_30s_.*\\.(tif|bil)$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  files[1]
}

extract_prism_wy <- function(sample_year, lats, lons) {
  wy_yyyymm    <- c(sprintf("%04d%02d", sample_year - 1, 10:12),
                    sprintf("%04d%02d", sample_year,     1:9))
  bio18_yyyymm <- sprintf("%04d%02d", sample_year, 6:8)

  bil_paths <- sapply(wy_yyyymm, find_prism_bil_r9)
  missing   <- wy_yyyymm[sapply(bil_paths, is.null)]
  if (length(missing) > 0)
    stop("Missing PRISM raster(s) for water year ", sample_year, ": ",
         paste(missing, collapse = ", "))

  stack        <- rast(unlist(bil_paths))
  names(stack) <- wy_yyyymm

  pts  <- vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326")
  pts  <- project(pts, crs(stack))
  vals <- as.matrix(extract(stack, pts)[, -1])

  bio12 <- rowSums(vals)
  bio15 <- apply(vals, 1, function(x) {
    if (any(is.na(x))) return(NA_real_)
    (sd(x) / mean(x)) * 100
  })
  bio18 <- rowSums(vals[, bio18_yyyymm, drop = FALSE])

  list(BIO12 = bio12, BIO15 = bio15, BIO18 = bio18)
}

# --- load unique coordinates ---

pheno        <- read.delim(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)
coords_range <- unique(pheno[, c("Locality", "Lat _dd.dddd", "Lon_dd.dddd")])
colnames(coords_range) <- c("Locality", "Lat", "Long")
coords_range <- coords_range[!is.na(coords_range$Lat) & !is.na(coords_range$Long), ]
coords_range <- coords_range[order(coords_range$Locality, coords_range$Lat, coords_range$Long), ]
rownames(coords_range) <- NULL
cat("\nSection 10: unique coordinates for range extraction:", nrow(coords_range), "\n")

# --- determine water-year range from available PRISM months ---
# Archive covers Oct 1926 – Sep 2024 → water years 1927–2024

wy_first <- 1927
wy_last  <- 2024
wy_years <- wy_first:wy_last
cat("Water years:", wy_first, "–", wy_last, "(", length(wy_years), "years )\n")

# --- pre-allocate output columns in the desired order ---

out <- coords_range
for (yr in wy_years)  out[[paste0("BIO12_", yr)]] <- NA_real_
for (yr in wy_years)  out[[paste0("BIO15_", yr)]] <- NA_real_
for (yr in wy_years)  out[[paste0("BIO18_", yr)]] <- NA_real_
fog_years <- fog_year_min:fog_year_max          # 2000:2025
for (yr in fog_years) out[[paste0("fog_",   yr)]] <- NA_real_

# --- extract PRISM ---

for (yr in wy_years) {
  cat("  Water year", yr, "...\n")
  tryCatch({
    pv <- extract_prism_wy(sample_year = yr,
                           lats        = coords_range$Lat,
                           lons        = coords_range$Long)
    out[[paste0("BIO12_", yr)]] <- pv$BIO12
    out[[paste0("BIO15_", yr)]] <- pv$BIO15
    out[[paste0("BIO18_", yr)]] <- pv$BIO18
  }, error = function(e) {
    warning("Could not extract PRISM data for water year ", yr, ": ", e$message)
  })
}

# --- extract fog (2000-2025, mean of Jun-Sep FLCC per year) ---

extract_fog_coords <- function(year, lats, lons) {
  paths <- file.path(fog_dir, sprintf("%d-%02d.tif", year, fog_months))
  names(paths) <- sprintf("%d-%02d", year, fog_months)
  missing_months <- names(paths)[!file.exists(paths)]
  if (length(missing_months) > 0)
    stop("Missing fog raster(s) for year ", year, ": ",
         paste(missing_months, collapse = ", "))
  pts <- vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326")
  monthly_vals <- matrix(
    sapply(unname(paths), function(p) {
      r        <- rast(p)
      pts_proj <- project(pts, crs(r))
      extract(r, pts_proj)[, 2]
    }),
    nrow = length(lats)
  )
  rowMeans(monthly_vals, na.rm = FALSE)
}

cat("Extracting fog values for range coordinates...\n")
for (yr in fog_years) {
  cat("  Fog year", yr, "...\n")
  tryCatch({
    out[[paste0("fog_", yr)]] <- extract_fog_coords(yr, coords_range$Lat, coords_range$Long)
  }, error = function(e) {
    warning("Could not extract fog data for year ", yr, ": ", e$message)
  })
}

# --- per-row averages across all available years ---

out$BIO12_mean <- rowMeans(out[, paste0("BIO12_", wy_years)],  na.rm = TRUE)
out$BIO15_mean <- rowMeans(out[, paste0("BIO15_", wy_years)],  na.rm = TRUE)
out$BIO18_mean <- rowMeans(out[, paste0("BIO18_", wy_years)],  na.rm = TRUE)
out$fog_mean   <- rowMeans(out[, paste0("fog_",   fog_years)], na.rm = TRUE)

# ============================================================
# 11. Kling annual fog from the real observed time series (1950-2014)
#     Kling/timeseries.tif : 65 layers, one per year, values are the
#     Jun-Sep mean fog frequency (same units/grid as the earlier
#     Kling_fog_normals.tif). This replaces the earlier trend+normal
#     reconstruction now that the actual annual data is available.
# ============================================================

cat("\nSection 11: Extracting real Kling fog annual values (1950-2014)...\n")

if (!file.exists(kling_ts_file)) {
  stop("Kling fog time series raster not found at ", kling_ts_file,
       ". Decompress Kling/timeseries.tif.zip first.")
}

r_kling_ts_r11  <- rast(kling_ts_file)
kling_years_r11 <- sort(as.integer(names(r_kling_ts_r11)))

pts_kling <- vect(data.frame(x = coords_range$Long, y = coords_range$Lat),
                  geom = c("x", "y"), crs = "EPSG:4326")
pts_kling <- project(pts_kling, crs(r_kling_ts_r11))

kling_vals <- as.matrix(
  extract(r_kling_ts_r11[[as.character(kling_years_r11)]], pts_kling)[, -1]
)

for (yr in kling_years_r11) out[[paste0("kling_fog_", yr)]] <- NA_real_
for (i in seq_along(kling_years_r11)) {
  out[[paste0("kling_fog_", kling_years_r11[i])]] <- kling_vals[, i]
}

out$kling_fog_mean <- rowMeans(out[, paste0("kling_fog_", kling_years_r11)], na.rm = TRUE)

n_kling <- sum(rowSums(is.na(out[, paste0("kling_fog_", kling_years_r11)])) == 0)
cat("Rows with complete kling_fog across 1950-2014:", n_kling, "/", nrow(out), "\n")

# --- save ---

# TEMPORARILY DISABLED — do not overwrite the existing output while
# developing/reviewing the _25yr environmental variable selection above.
# write.csv(out, range_output, row.names = FALSE)
cat("\nRange extraction complete. Saved to:", range_output, "\n")
cat("Dimensions:", nrow(out), "rows ×", ncol(out), "columns\n")

n_complete <- sum(rowSums(is.na(out[, paste0("BIO12_", wy_years)])) == 0)
cat("Rows with complete BIO12 across all years:", n_complete, "/", nrow(out), "\n")
n_fog <- sum(rowSums(is.na(out[, paste0("fog_", fog_years)])) == 0)
cat("Rows with complete fog across 2000-2025:", n_fog, "/", nrow(out), "\n")

} # end if (RUN_FULL_EXTRACTION)

# ============================================================
# 12. BCM climatic water deficit (CWD_25yr)
#
#     Adds one new column to emma_gee_results_final.csv: CWD_25yr, the
#     fixed-reference-window mean annual CWD (mm/water-year), computed
#     the same way as BIO12_25yr/BIO15_25yr/BIO18_25yr above — averaged
#     once per unique sample coordinate, then joined back onto every
#     sample sharing that coordinate.
#
#     Source: BCMv8 statewide monthly CWD rasters (same data used by
#     env/extract_BCM_CWD.R), NAD83 California Teale Albers (EPSG:3310),
#     270 m, at /private/groups/meyerlab/esteiger/env/BCM/.
#
#     Window caveat (resolved 2026-10-22): BIO12/15/18_25yr average water
#     years 2000-2024 (PRISM archive runs through Sep 2024). The BCM CWD
#     archive on disk previously only had monthly rasters through WY2022;
#     cwd_WY2023.zip, cwd_WY2024.zip, and cwd_WY2025.zip have since been
#     downloaded from ScienceBase (same BCMv8 statewide CWD collection,
#     DOI 10.5066/P9PT36UI — see env/extract_BCM_CWD.R for item IDs), so
#     CWD_25yr is now a 26-water-year mean (2000-2025), matching/exceeding
#     the other "_25yr" columns' coverage rather than trailing them.
# ============================================================

cat("\n=== Section 12: BCM CWD_25yr extraction ===\n\n")

# Unzipping+summing 26 water years of BCM rasters takes a while (~20 min
# for the original 23 years). Once CWD_25yr is in output_file for every
# sample, skip re-doing that on every re-run — flip to TRUE (or delete
# the CWD_25yr column from the csv) to force a re-extraction.
cwd_bcm_dir   <- "/private/groups/meyerlab/esteiger/env/BCM"
cwd_ref_years <- 2000:2025   # see window caveat above

# --- helpers (same logic as env/extract_BCM_CWD.R) — defined unconditionally
#     (cheap) so Section 14 can reuse them for a handful of new coordinates
#     even when the full extraction below is skipped ---

wy_files_12 <- function(yr) {
  mon_names <- c("oct","nov","dec","jan","feb","mar","apr","may","jun","jul","aug","sep")
  years     <- c(rep(yr - 1, 3), rep(yr, 9))
  sprintf("cwd%d%s.asc", years, mon_names)
}

get_wy_zip_12 <- function(yr) {
  if      (yr <= 2009) file.path(cwd_bcm_dir, "cwd_WY2000_09.zip")
  else if (yr <= 2020) file.path(cwd_bcm_dir, "cwd_WY2010_20.zip")
  else if (yr <= 2025) file.path(cwd_bcm_dir, sprintf("cwd_WY%d.zip", yr))
  else stop(sprintf("No BCM CWD zip available for WY%d", yr))
}

# terra 1.5.21 buffer extraction does not search outward from a NoData
# cell, so NA points get a focal-mean fill (grows the search window)
# before falling back to leaving them NA.
fill_na_focal_12 <- function(raster, pts_3310, vals, windows = c(3, 7, 15)) {
  na_idx <- which(is.na(vals))
  if (length(na_idx) == 0) return(vals)
  na_ext <- ext(pts_3310[na_idx, ]) + 5000
  r_crop <- crop(raster, na_ext)
  for (w in windows) {
    still_na <- na_idx[is.na(vals[na_idx])]
    if (length(still_na) == 0) break
    r_filled <- focal(r_crop, w = matrix(1, w, w), fun = mean, na.rm = TRUE, na.policy = "only")
    vals[still_na] <- extract(r_filled, pts_3310[still_na, ])[, 2]
  }
  vals
}

sum_wy_cwd_12 <- function(zip_path, file_names, study_ext) {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE))
  cat(sprintf("  Unzipping %d files from %s...\n", length(file_names), basename(zip_path)))
  unzip(zip_path, files = file_names, exdir = tmpdir)
  rasters <- lapply(file.path(tmpdir, file_names), function(f) {
    r <- rast(f)
    crs(r) <- "EPSG:3310"
    crop(r, study_ext)
  })
  Reduce("+", rasters)
}

# Shared by both Section 12 (all sample coordinates) and Section 14 (a
# handful of new ones): mean CWD over cwd_ref_years at arbitrary Lat/Long.
extract_cwd_25yr <- function(lats, longs) {
  xy <- data.frame(Lat = lats, Long = longs)
  pts_wgs84 <- vect(xy, geom = c("Long", "Lat"), crs = "EPSG:4326")
  pts_3310  <- project(pts_wgs84, "EPSG:3310")
  study_ext <- ext(pts_3310) + 5000

  mat <- matrix(NA_real_, nrow = nrow(xy), ncol = length(cwd_ref_years),
                dimnames = list(NULL, cwd_ref_years))
  for (yr in cwd_ref_years) {
    zip   <- get_wy_zip_12(yr)
    files <- wy_files_12(yr)
    cat(sprintf("WY%d: reading from %s\n", yr, basename(zip)))
    zip_contents <- unzip(zip, list = TRUE)$Name
    missing <- files[!files %in% zip_contents]
    if (length(missing) > 0) {
      warning(sprintf("WY%d: %d file(s) missing from zip, skipping: %s",
                      yr, length(missing), paste(missing, collapse = ", ")))
      next
    }
    cwd_wy <- sum_wy_cwd_12(zip, files, study_ext)
    vals   <- extract(cwd_wy, pts_3310)[, 2]
    vals   <- fill_na_focal_12(cwd_wy, pts_3310, vals)
    mat[, as.character(yr)] <- vals
    cat(sprintf("  WY%d done: mean = %.1f mm\n", yr, mean(vals, na.rm = TRUE)))
  }
  rowMeans(mat, na.rm = TRUE)
}

RUN_CWD_EXTRACTION <- {
  existing_check <- tryCatch(read.csv(output_file), error = function(e) NULL)
  is.null(existing_check) || !("CWD_25yr" %in% names(existing_check)) ||
    anyNA(existing_check$CWD_25yr)
}

if (RUN_CWD_EXTRACTION) {

coords12 <- read.csv(coords_file)
cat("Loaded", nrow(coords12), "sample coordinates\n")

uniq_xy12 <- unique(coords12[, c("Lat", "Long")])
rownames(uniq_xy12) <- NULL
cat("Unique coordinates:", nrow(uniq_xy12), "\n")

# --- extract annual CWD at the unique coordinates, one water year at a time ---

uniq_xy12$CWD_25yr <- extract_cwd_25yr(uniq_xy12$Lat, uniq_xy12$Long)
cat(sprintf("\nCWD_25yr (WY%d-%d mean): %d / %d coordinates non-NA\n",
            min(cwd_ref_years), max(cwd_ref_years),
            sum(!is.na(uniq_xy12$CWD_25yr)), nrow(uniq_xy12)))

# --- join onto every sample sharing that coordinate, merge into the
#     existing emma_gee_results_final.csv, and write it back ---
#
#     CWD_25yr depends only on (Lat, Long), never on SampleID, so it's
#     assigned via a direct keyed lookup against the ORIGINAL row order
#     -- never reconstructing/reordering the dataframe through merge() +
#     match(existing$SampleID, ...). That reorder-by-SampleID approach
#     silently corrupted OTHER columns (fog_modern, BIO18_modern, etc.)
#     for any SampleID that appears more than once in the file (found
#     2026-10-22: SampleIDs 6083 and 8071 each have two rows sharing one
#     coordinate; match() mapped both of a duplicated ID's rows onto
#     whichever one came first in merged's internal order, discarding
#     the other row's real values). existing$SampleID is NOT guaranteed
#     unique (confirmed duplicates exist), so never key a reorder on it.

existing <- read.csv(output_file)
existing$CWD_25yr <- NULL   # drop any stale copy from a previous run
key_existing <- paste(round(existing$Lat, 8), round(existing$Long, 8))
key_uniq     <- paste(round(uniq_xy12$Lat, 8), round(uniq_xy12$Long, 8))
existing$CWD_25yr <- uniq_xy12$CWD_25yr[match(key_existing, key_uniq)]
merged <- existing

write.csv(merged, output_file, row.names = FALSE)
cat("\nCWD_25yr added to", output_file, "\n")
cat("  Non-NA CWD_25yr:", sum(!is.na(merged$CWD_25yr)), "/", nrow(merged), "samples\n")

} else {
  cat("CWD_25yr already present (no NAs) in", output_file, "— skipping re-extraction.\n",
      "Set RUN_CWD_EXTRACTION <- TRUE above to force a re-run.\n")
}

# ============================================================
# 13. Extract env vars for observations previously missing from the
#     coordinate-extraction pipeline
#
#     Two sources of "new" modern (year>2000) observations not yet in
#     emma_gee_results_final.csv, found while reconciling with the
#     undergrad's ray-floret-count protocol (2026-09-05):
#       (a) gtinatandvars.csv (the undergrad's original phenotype file,
#           received 2026-09-06 — before that this section used
#           updatedvarsandpheno.csv, a since-superseded bridging copy of
#           the same 224 specimens) rows with a valid recordNumber_
#           CollectorNumber + Lat/Long that simply were never added to
#           all_coord_gvtp.csv (e.g. specimens 8503-8505, 17974) — real
#           specimens, just missing from the coordinate-tracking file,
#           not "weird" exclusions as their absence from the models had
#           suggested. Rows with a date stored in place of a collector
#           number (e.g. "2001-07-10") are skipped — not a real ID.
#       (b) inatrays.csv rows (iNat-only observations, no physical
#           specimen, ray-floret count only) with a usable iNat_ray_count.
#     For each new coordinate, first tries a match (~50m tolerance)
#     against already-extracted coordinates in emma_gee_results_final.csv
#     — cheap, no raster work; most (43/46 as of 2026-09-05) resolve this
#     way since the same field sites were often visited/photographed more
#     than once. Only genuinely new coordinates (3/46) get a fresh
#     extraction, via the same helpers as Sections 5-7/12 but redefined
#     minimally here so this section doesn't depend on RUN_FULL_EXTRACTION.
#     Self-limiting: an already-processed SampleID is excluded from
#     "new_samples", so re-running this script after the first successful
#     pass finds nothing left to do.
# ============================================================

cat("\n=== Section 13: env extraction for previously-uncoordinated observations ===\n\n")

pheno_updated_file <- "/private/groups/meyerlab/esteiger/pheno/gtinatandvars.csv"
inatrays_file      <- "/private/groups/meyerlab/esteiger/pheno/inatrays.csv"

find_prism_bil_13 <- function(yyyymm) {
  dir_path <- file.path(prism_dir, yyyymm)
  files <- list.files(dir_path, pattern = "_30s_.*\\.(tif|bil)$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  files[1]
}

extract_fog_seasonal_13 <- function(year, lats, lons) {
  paths <- file.path(fog_dir, sprintf("%d-%02d.tif", year, fog_months))
  names(paths) <- sprintf("%d-%02d", year, fog_months)
  missing_months <- names(paths)[!file.exists(paths)]
  if (length(missing_months) > 0) {
    warning("Missing fog raster(s) for year ", year, ": ", paste(missing_months, collapse = ", "))
    return(rep(NA_real_, length(lats)))
  }
  pts <- vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326")
  monthly_vals <- matrix(
    sapply(unname(paths), function(p) {
      r <- rast(p); pts_proj <- project(pts, crs(r)); extract(r, pts_proj)[, 2]
    }),
    nrow = length(lats)
  )
  rowMeans(monthly_vals, na.rm = FALSE)
}

extract_fog_25yr_13 <- function(lats, lons) {
  mat <- matrix(NA_real_, nrow = length(lats), ncol = length(ref_fog_years))
  for (i in seq_along(ref_fog_years)) mat[, i] <- extract_fog_seasonal_13(ref_fog_years[i], lats, lons)
  rowMeans(mat, na.rm = TRUE)
}

extract_prism_year_13 <- function(sample_year, lats, lons) {
  wy_yyyymm    <- c(sprintf("%04d%02d", sample_year - 1, 10:12), sprintf("%04d%02d", sample_year, 1:9))
  bio18_yyyymm <- sprintf("%04d%02d", sample_year, 6:8)
  bil_paths <- sapply(wy_yyyymm, find_prism_bil_13)
  missing   <- wy_yyyymm[sapply(bil_paths, is.null)]
  if (length(missing) > 0) {
    warning("Missing PRISM raster(s) for water year ", sample_year, ": ", paste(missing, collapse = ", "))
    n <- length(lats)
    return(list(BIO12 = rep(NA_real_, n), BIO15 = rep(NA_real_, n), BIO18 = rep(NA_real_, n)))
  }
  stack        <- rast(unlist(bil_paths))
  names(stack) <- wy_yyyymm
  pts <- vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326")
  pts <- project(pts, crs(stack))
  vals <- as.matrix(extract(stack, pts)[, -1])
  bio12 <- rowSums(vals)
  bio15 <- apply(vals, 1, function(x) if (any(is.na(x))) NA_real_ else (sd(x) / mean(x)) * 100)
  bio18 <- rowSums(vals[, bio18_yyyymm, drop = FALSE])
  list(BIO12 = bio12, BIO15 = bio15, BIO18 = bio18)
}

extract_bio_25yr_13 <- function(lats, lons) {
  n <- length(lats)
  bio12 <- bio15 <- bio18 <- matrix(NA_real_, nrow = n, ncol = length(ref_bio_years))
  for (i in seq_along(ref_bio_years)) {
    pv <- extract_prism_year_13(ref_bio_years[i], lats, lons)
    bio12[, i] <- pv$BIO12; bio15[, i] <- pv$BIO15; bio18[, i] <- pv$BIO18
  }
  list(BIO12_25yr = rowMeans(bio12, na.rm = TRUE), BIO15_25yr = rowMeans(bio15, na.rm = TRUE),
       BIO18_25yr = rowMeans(bio18, na.rm = TRUE))
}

existing13 <- read.csv(output_file, stringsAsFactors = FALSE)
existing13$SampleID <- trimws(as.character(existing13$SampleID))

# --- (a) previously-unmatched gtinatandvars.csv rows ---
pv13 <- read.csv(pheno_updated_file, stringsAsFactors = FALSE)
pv13$recordNumber_CollectorNumber <- trimws(as.character(pv13$recordNumber_CollectorNumber))
pv13$year <- suppressWarnings(as.numeric(pv13$year))
pv13_new <- pv13[!is.na(pv13$year) & pv13$year > 2000 &
                 !is.na(pv13$Lat._dd.dddd) & !is.na(pv13$Lon_dd.dddd) &
                 !grepl("-", pv13$recordNumber_CollectorNumber) &
                 !(pv13$recordNumber_CollectorNumber %in% existing13$SampleID), ]
new_a <- data.frame(SampleID = pv13_new$recordNumber_CollectorNumber,
                    Lat = pv13_new$Lat._dd.dddd, Long = pv13_new$Lon_dd.dddd,
                    Year = pv13_new$year)

# --- (b) inatrays.csv rows with a usable ray count ---
inat13 <- read.csv(inatrays_file, stringsAsFactors = FALSE)
inat13$inat_Specimen <- trimws(as.character(inat13$inat_Specimen))
inat13_new <- inat13[!is.na(inat13$iNat_ray_count) & !is.na(inat13$Year) & inat13$Year > 2000 &
                      !is.na(inat13$latitude) & !is.na(inat13$longitude) &
                      !(inat13$inat_Specimen %in% existing13$SampleID), ]
new_b <- data.frame(SampleID = inat13_new$inat_Specimen,
                    Lat = inat13_new$latitude, Long = inat13_new$longitude,
                    Year = inat13_new$Year)

new_samples <- unique(rbind(new_a, new_b))
new_samples <- new_samples[!duplicated(new_samples$SampleID), ]
cat("New samples needing env data:", nrow(new_samples), "\n")

if (nrow(new_samples) > 0) {
  uniq_new_xy <- unique(new_samples[, c("Lat", "Long")])
  rownames(uniq_new_xy) <- NULL
  cat("Unique new coordinates:", nrow(uniq_new_xy), "\n")

  match_idx <- sapply(seq_len(nrow(uniq_new_xy)), function(i) {
    hit <- which(abs(existing13$Lat - uniq_new_xy$Lat[i]) < 0.0005 &
                 abs(existing13$Long - uniq_new_xy$Long[i]) < 0.0005)
    if (length(hit) > 0) hit[1] else NA_integer_
  })

  uniq_new_xy$fog_25yr   <- NA_real_
  uniq_new_xy$BIO18_25yr <- NA_real_
  uniq_new_xy$BIO15_25yr <- NA_real_
  uniq_new_xy$CWD_25yr   <- NA_real_

  has_match <- !is.na(match_idx)
  cat("Recovered via existing coordinate match:", sum(has_match), "\n")
  uniq_new_xy$fog_25yr[has_match]   <- existing13$fog_25yr[match_idx[has_match]]
  uniq_new_xy$BIO18_25yr[has_match] <- existing13$BIO18_25yr[match_idx[has_match]]
  uniq_new_xy$BIO15_25yr[has_match] <- existing13$BIO15_25yr[match_idx[has_match]]
  uniq_new_xy$CWD_25yr[has_match]   <- existing13$CWD_25yr[match_idx[has_match]]

  needs_fresh <- !has_match
  cat("Needing fresh raster extraction:", sum(needs_fresh), "\n")

  if (any(needs_fresh)) {
    fl <- uniq_new_xy$Lat[needs_fresh]; fg <- uniq_new_xy$Long[needs_fresh]

    cat("  fog_25yr...\n")
    uniq_new_xy$fog_25yr[needs_fresh] <- extract_fog_25yr_13(fl, fg)

    cat("  BIO15/18_25yr...\n")
    bio_fresh <- extract_bio_25yr_13(fl, fg)
    uniq_new_xy$BIO15_25yr[needs_fresh] <- bio_fresh$BIO15_25yr
    uniq_new_xy$BIO18_25yr[needs_fresh] <- bio_fresh$BIO18_25yr

    cat("  CWD_25yr...\n")
    uniq_new_xy$CWD_25yr[needs_fresh] <- extract_cwd_25yr(fl, fg)
  }

  new_samples_env <- merge(new_samples, uniq_new_xy, by = c("Lat", "Long"), all.x = TRUE)
  for (col in setdiff(names(existing13), names(new_samples_env))) new_samples_env[[col]] <- NA
  new_samples_env <- new_samples_env[, names(existing13)]

  updated13 <- rbind(existing13, new_samples_env)
  write.csv(updated13, output_file, row.names = FALSE)
  cat("\nAppended", nrow(new_samples_env), "new samples to", output_file, "\n")
  cat("New total rows:", nrow(updated13), "\n")
} else {
  cat("No new samples to extract — nothing to do.\n")
}

# ============================================================
# 13b. Backfill _germyr values for any row still missing them
#
#     Independent of the "new samples" logic above (which only fires the
#     first time a SampleID appears) — samples this section added in a
#     PRIOR run already have fog_25yr/BIO18_25yr/BIO15_25yr/CWD_25yr but
#     were never given _germyr values, since they have no row at all in
#     the original germyr extraction (Section 6), having not existed in
#     all_coord_gvtp.csv when that ran. Discovered 2026-09-06 via
#     modern_data.R's turbine-effect model (which needs germyr, not
#     _25yr, covariates): 31 of the ~35 relevant specimens at the 3
#     turbine-analysis localities had NA fog/BIO18/BIO15_germyr for
#     exactly this reason, collapsing the "before" sample from ~35 to 4.
#     Re-derives each affected row's OWN year from the specimen files
#     (gtinatandvars.csv or inatrays.csv) to extract at the right time,
#     reusing the per-year extraction helpers already defined above.
#     Self-limiting/idempotent: only rows with NA fog_germyr are touched.
# ============================================================

cat("\n=== Section 13b: backfill _germyr for rows still missing it ===\n\n")

existing13b <- read.csv(output_file, stringsAsFactors = FALSE)
existing13b$SampleID <- trimws(as.character(existing13b$SampleID))
gap_idx <- which(is.na(existing13b$fog_germyr))
cat("Rows missing fog_germyr:", length(gap_idx), "\n")

if (length(gap_idx) > 0) {
  # recover each gap row's own collection year from whichever specimen
  # file has it (gtinatandvars.csv for physical/pheno_file specimens,
  # inatrays.csv for iNat-only ones)
  pv13b <- read.csv(pheno_updated_file, stringsAsFactors = FALSE)
  pv13b$recordNumber_CollectorNumber <- trimws(as.character(pv13b$recordNumber_CollectorNumber))
  year_by_id <- setNames(suppressWarnings(as.numeric(pv13b$year)), pv13b$recordNumber_CollectorNumber)

  inat13b <- read.csv(inatrays_file, stringsAsFactors = FALSE)
  inat13b$inat_Specimen <- trimws(as.character(inat13b$inat_Specimen))
  year_by_id_inat <- setNames(inat13b$Year, inat13b$inat_Specimen)

  gap_ids <- existing13b$SampleID[gap_idx]
  gap_year <- ifelse(!is.na(year_by_id[gap_ids]), year_by_id[gap_ids], year_by_id_inat[gap_ids])
  existing13b$Year[gap_idx] <- ifelse(is.na(existing13b$Year[gap_idx]), gap_year, existing13b$Year[gap_idx])

  still_missing_year <- gap_idx[is.na(existing13b$Year[gap_idx])]
  if (length(still_missing_year) > 0) {
    cat("  (", length(still_missing_year), "row(s) have no recoverable year — left NA)\n")
  }
  gap_idx <- setdiff(gap_idx, still_missing_year)

  for (yr in sort(unique(existing13b$Year[gap_idx]))) {
    idx <- gap_idx[existing13b$Year[gap_idx] == yr]
    lats <- existing13b$Lat[idx]; longs <- existing13b$Long[idx]
    cat(sprintf("  germyr for Year %d (%d sample(s))...\n", yr, length(idx)))
    if (yr >= fog_year_min && yr <= fog_year_max) {
      existing13b$fog_germyr[idx] <- extract_fog_seasonal_13(yr, lats, longs)
    }
    pv <- extract_prism_year_13(yr, lats, longs)
    existing13b$BIO12_germyr[idx] <- pv$BIO12
    existing13b$BIO15_germyr[idx] <- pv$BIO15
    existing13b$BIO18_germyr[idx] <- pv$BIO18
  }
  write.csv(existing13b, output_file, row.names = FALSE)
  cat(sprintf("\ngermyr backfill: %d / %d previously-missing rows now non-NA fog_germyr\n",
              sum(!is.na(existing13b$fog_germyr[gap_idx])), length(gap_idx)))
  cat("Written to", output_file, "\n")
} else {
  cat("No rows missing _germyr — nothing to do.\n")
}

# ============================================================
# 14. Collinearity check: fog/BIO18/BIO15/CWD_25yr (modern_data.R's vars)
#
#     Locality-capped (15/locality, same method + seed as Section 5b.3),
#     coordinate-deduplicated Spearman correlations among the four _25yr
#     variables modern_data.R actually models. Self-contained — reads the
#     phenotype tsv fresh and pulls already-extracted values back out of
#     emma_gee_results_final.csv, so it doesn't require re-running the
#     expensive bioclim-analog extraction in Section 5b just to add CWD
#     to that analysis.
#
#     Result (2026-09-05, n=63 unique locality-capped coordinates): no
#     pair exceeds the 0.7 threshold used elsewhere in this script.
#     Strongest pair: BIO18_25yr x CWD_25yr, rho = -0.60 (mechanistically
#     expected — CWD is water deficit, BIO18 is summer precip).
#
#     A naive check with no locality cap and no coordinate dedup (all
#     1275 GEN+PHENO sample rows) had instead suggested fog x BIO18 was
#     collinear (rho = -0.68 to -0.76 depending on which sample). That
#     was an artifact of heavily oversampled localities (e.g. Strauss
#     n=50, VAFB_Lion_head n=48 in the phenotype set) dominating the
#     correlation — it disappears once each locality is capped and
#     coordinates are deduplicated, confirming the same 15/locality-capped
#     methodology used elsewhere in this script.
# ============================================================

cat("\n=== Section 13: fog/BIO18/BIO15/CWD_25yr collinearity check ===\n\n")

MAX_PER_LOCALITY_CC     <- 15    # matches Section 5b.3
COLLINEAR_THRESH_25yr_CC <- 0.7  # matches Section 5b.4
modeling_vars <- c("fog_25yr", "BIO18_25yr", "BIO15_25yr", "CWD_25yr")

pheno_cc <- read.delim(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)
pheno_cc <- pheno_cc[, c("Locality", "Lat _dd.dddd", "Lon_dd.dddd")]
colnames(pheno_cc) <- c("Locality", "Lat", "Long")
pheno_cc <- pheno_cc[!is.na(pheno_cc$Locality) & pheno_cc$Locality != "" &
                     !is.na(pheno_cc$Lat) & !is.na(pheno_cc$Long), ]

set.seed(42)
capped_cc <- do.call(rbind, lapply(split(pheno_cc, pheno_cc$Locality), function(d) {
  if (nrow(d) > MAX_PER_LOCALITY_CC) d[sample(nrow(d), MAX_PER_LOCALITY_CC), ] else d
}))
rownames(capped_cc) <- NULL

uniq_xy_cc <- unique(capped_cc[, c("Lat", "Long")])
rownames(uniq_xy_cc) <- NULL
cat("Locality-capped, deduplicated coordinates:", nrow(uniq_xy_cc),
    "(from", nrow(capped_cc), "capped rows,", nrow(pheno_cc), "total phenotype rows)\n")

env_cc      <- read.csv(output_file)
env_uniq_cc <- unique(env_cc[, c("Lat", "Long", modeling_vars)])
merged_cc   <- merge(uniq_xy_cc, env_uniq_cc, by = c("Lat", "Long"), all.x = TRUE)
cat("Matched env data for", sum(complete.cases(merged_cc[, modeling_vars])), "/",
    nrow(merged_cc), "unique coordinates\n\n")

cor_mat_cc <- cor(merged_cc[, modeling_vars], method = "spearman", use = "pairwise.complete.obs")
cat("Spearman correlations (locality-capped, n =", nrow(merged_cc), "unique coords):\n")
print(round(cor_mat_cc, 2))

flagged_cc <- which(abs(cor_mat_cc) >= COLLINEAR_THRESH_25yr_CC & upper.tri(cor_mat_cc), arr.ind = TRUE)
if (nrow(flagged_cc) == 0) {
  cat("\nNo pair exceeds the |rho| >=", COLLINEAR_THRESH_25yr_CC, "threshold.\n")
} else {
  cat("\nFlagged pairs (|rho| >=", COLLINEAR_THRESH_25yr_CC, "):\n")
  for (k in seq_len(nrow(flagged_cc))) {
    cat(" ", rownames(cor_mat_cc)[flagged_cc[k, 1]], "x", colnames(cor_mat_cc)[flagged_cc[k, 2]],
        "=", round(cor_mat_cc[flagged_cc[k, 1], flagged_cc[k, 2]], 2), "\n")
  }
}

cor_df_cc <- data.frame(
  var1 = rep(rownames(cor_mat_cc), ncol(cor_mat_cc)),
  var2 = rep(colnames(cor_mat_cc), each = nrow(cor_mat_cc)),
  rho  = as.vector(cor_mat_cc)
)

jpeg(file.path(plots_dir, "env_selection_correlation_heatmap_25yr_modeling_vars.jpg"),
     width = 800, height = 700)
print(
  ggplot(cor_df_cc, aes(x = var2, y = var1, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 5) +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
    theme_bw(base_size = 14) +
    theme(axis.title = element_blank()) +
    ggtitle(sprintf("modern_data.R variables: locality-capped (n=%d) Spearman correlations",
                    nrow(merged_cc)))
)
dev.off()
cat("\nHeatmap written to",
    file.path(plots_dir, "env_selection_correlation_heatmap_25yr_modeling_vars.jpg"), "\n")

# ============================================================
# 15. Rolling 20-year (inclusive of germination year) env variables:
#     _modern (fog/BIO18/BIO15, for modern_data.R) and _historic
#     (Kling/BIO18/BIO15/CWD, for historic_data.R)
#
#     Replaces the earlier FIXED-reference-window _25yr/_65yr approach
#     (Sections 5b/7/12/13 above, no longer used by modern_data.R/
#     historic_data.R, but left in place rather than deleted) with a
#     window tied to each specimen's OWN germination year:
#     [germYear - (ROLL_WINDOW-1), germYear], ROLL_WINDOW=20 years,
#     INCLUDING the germination year itself (window length matched to
#     modern_data.R's own CWD window, user decision 2026-10-26). Not every
#     one of those 20 years needs to be present -- see MIN_WINDOW_COVERAGE
#     below -- so a specimen isn't dropped just for 1-2 missing years.
#
#     CWD is now included (as of 2026-10-26), via
#     env/crop_cache_bcm_cwd_annual.R's pre-cropped per-water-year raster
#     cache (env/BCM/cwd_annual_cropped/) rather than re-unzipping a
#     statewide BCM zip on every germYear window -- that per-window unzip
#     cost is what got CWD dropped from this section originally, when the
#     cache only covered WY2000-2022. The cache was later extended back to
#     WY1900 (see that script), comfortably covering the full historic
#     specimen record (germYear 1927-2024), so CWD is no longer the
#     binding constraint for historic_data.R the way Kling's narrower
#     1950-2014 range still is.
#
#     Coverage constraints (via safe_rowMeans()'s MIN_WINDOW_COVERAGE, not
#     an all-or-nothing requirement that every one of the ROLL_WINDOW
#     years be present): a specimen's germYear must fall within
#     [source_min + ROLL_WINDOW - 1, source_max] (+/- GATE_BUFFER) for
#     that source:
#       fog_modern:      fog (Werner): 2000-2025
#       BIO*_modern/historic: PRISM: 1927-2024 -- wide enough to never
#                         bind for either modern_data.R or historic_data.R
#       Kling_historic:  Kling: 1950-2014 -- the binding constraint for
#                         historic_data.R's fog/Kling variable
#       CWD_historic:    BCM cache: see cwd_cache_years above (adapts
#                         automatically to however much of WY1900-2025 has
#                         been downloaded/cached)
# ============================================================

cat("\n=== Section 15: rolling 10yr (incl. germ year) env variables ===\n\n")

ROLL_WINDOW <- 20  # years, inclusive of the germination year itself

r_kling_ts_15  <- rast(kling_ts_file)
kling_years_15 <- sort(as.integer(names(r_kling_ts_15)))

extract_kling_seasonal_15 <- function(year, lats, lons) {
  yr_chr <- as.character(year)
  if (!(yr_chr %in% names(r_kling_ts_15))) return(rep(NA_real_, length(lats)))
  pts <- vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326")
  pts <- project(pts, crs(r_kling_ts_15))
  extract(r_kling_ts_15[[yr_chr]], pts)[, 2]
}

# --- CWD, now feasible here thanks to env/crop_cache_bcm_cwd_annual.R's
#     pre-cropped per-water-year cache (env/BCM/cwd_annual_cropped/
#     cwd_WY<year>.tif, one small raster per year rather than a statewide
#     zip needing a fresh unzip+crop+sum every time a germYear window
#     touches it -- this is exactly the cost that got CWD dropped from this
#     section in the first place, per the comment above). Cache spans
#     whatever decade zips have been downloaded/processed; computed
#     dynamically below so this adapts automatically as more years are added.
cwd_cache_dir   <- "/private/groups/meyerlab/esteiger/env/BCM/cwd_annual_cropped"
cwd_cache_files <- list.files(cwd_cache_dir, pattern = "^cwd_WY[0-9]{4}\\.tif$")
cwd_cache_years <- sort(as.integer(gsub("[^0-9]", "", cwd_cache_files)))
cat(sprintf("CWD annual cache: %d year(s) available (%s)\n", length(cwd_cache_years),
            if (length(cwd_cache_years) > 0) paste(range(cwd_cache_years), collapse = "-") else "none"))

extract_cwd_annual_15 <- function(year, lats, lons) {
  if (!(year %in% cwd_cache_years)) return(rep(NA_real_, length(lats)))
  r <- rast(file.path(cwd_cache_dir, sprintf("cwd_WY%d.tif", year)))
  pts <- vect(data.frame(x = lons, y = lats), geom = c("x", "y"), crs = "EPSG:4326")
  pts <- project(pts, crs(r))
  extract(r, pts)[, 2]
}

# --- assemble the full specimen list (phenotype file + iNat-only rows,
#     same sourcing/dedup logic as Section 13) with each one's own
#     germination year and coordinates ---

pheno15 <- read.csv(pheno_updated_file, stringsAsFactors = FALSE)
pheno15$recordNumber_CollectorNumber <- trimws(as.character(pheno15$recordNumber_CollectorNumber))
pheno15$year <- suppressWarnings(as.numeric(pheno15$year))
pheno15_ok <- pheno15[!is.na(pheno15$year) & !is.na(pheno15$Lat._dd.dddd) & !is.na(pheno15$Lon_dd.dddd) &
                      !grepl("-", pheno15$recordNumber_CollectorNumber), ]
samples15 <- data.frame(SampleID = pheno15_ok$recordNumber_CollectorNumber,
                        Lat = pheno15_ok$Lat._dd.dddd, Long = pheno15_ok$Lon_dd.dddd,
                        GermYear = pheno15_ok$year)

inat15 <- read.csv(inatrays_file, stringsAsFactors = FALSE)
inat15$inat_Specimen <- trimws(as.character(inat15$inat_Specimen))
inat15_ok <- inat15[!is.na(inat15$iNat_ray_count) & !is.na(inat15$Year) &
                    !is.na(inat15$latitude) & !is.na(inat15$longitude) &
                    !(inat15$inat_Specimen %in% samples15$SampleID), ]
samples15 <- rbind(samples15, data.frame(SampleID = inat15_ok$inat_Specimen,
                                          Lat = inat15_ok$latitude, Long = inat15_ok$longitude,
                                          GermYear = inat15_ok$Year))
samples15 <- samples15[!duplicated(samples15$SampleID), ]
cat("Total specimens (any year, valid coords):", nrow(samples15), "\n")

# --- compute rolling-window means, one unique GermYear at a time (each
#     GermYear's window is fetched fresh even where it overlaps a
#     neighboring GermYear's window — simpler than caching per-year
#     rasters across groups, and fog/PRISM reads are fast enough that
#     this isn't a real cost, unlike the BCM CWD extraction) ---

# fog/Kling and BIO12/15/18 come from different sources with different
# coverage (fog: 2000-2025; Kling: 1950-2014; PRISM/BIO: 1927-2024) — a
# specimen whose window fits one but not the other still gets that one
# variable computed (each variable's own NA gate is independent). Found
# 2026-09-06: a single shared gate had been needlessly capping
# BIO18_historic/BIO15_historic (PRISM-based, wide coverage) to Kling's
# much narrower range, undercutting a mediation test that only needed
# BIO18, not Fog/Kling.
# A specimen's window can have a handful of missing years even when its
# GermYear itself is in-range (found 2026-09-06: 4 of 16 Kling-eligible
# historic specimens had a sporadic gap inside an otherwise-valid 1950-2014
# window). Require at least MIN_WINDOW_COVERAGE of the window's years
# rather than every single one, so a specimen isn't dropped over 1-2
# missing years out of 20.
MIN_WINDOW_COVERAGE <- 0.9
safe_rowMeans <- function(mat, min_frac = MIN_WINDOW_COVERAGE) {
  means <- rowMeans(mat, na.rm = TRUE)
  n_valid <- rowSums(!is.na(mat))
  means[n_valid < ceiling(min_frac * ncol(mat))] <- NA
  means
}
# The eligibility gate (below) must allow the same slack safe_rowMeans
# tolerates within a window, or a specimen whose window trails just past
# the archive's last available year never even reaches safe_rowMeans to
# be judged on its actual coverage. Found 2026-09-07: a germYear=2026
# specimen's window (2007-2026) is missing only 2026 itself — 19/20 =
# 95% coverage, comfortably above MIN_WINDOW_COVERAGE — but was being
# hard-excluded by a gate that only checked GermYear <= archive's last
# year, with no tolerance at all. GATE_BUFFER lets the gate itself admit
# a specimen whose window could still clear MIN_WINDOW_COVERAGE even if
# every one of those extra years turns out missing; safe_rowMeans then
# makes the real per-specimen call from actual per-year data.
GATE_BUFFER <- ROLL_WINDOW - ceiling(MIN_WINDOW_COVERAGE * ROLL_WINDOW)  # 20*0.9 -> 2 years

compute_rolling_15 <- function(samples, source, fog_min, fog_max, bio_min, bio_max, cwd_min, cwd_max) {
  n <- nrow(samples)
  val_fog <- val_bio12 <- val_bio15 <- val_bio18 <- val_cwd <- rep(NA_real_, n)
  ok_fog <- !is.na(samples$GermYear) & samples$GermYear >= fog_min - GATE_BUFFER & samples$GermYear <= fog_max + GATE_BUFFER
  ok_bio <- !is.na(samples$GermYear) & samples$GermYear >= bio_min - GATE_BUFFER & samples$GermYear <= bio_max + GATE_BUFFER
  ok_cwd <- !is.na(samples$GermYear) & samples$GermYear >= cwd_min - GATE_BUFFER & samples$GermYear <= cwd_max + GATE_BUFFER
  ok_any <- ok_fog | ok_bio | ok_cwd
  cat(sprintf("  %s: fog/Kling in-range=%d, BIO12/15/18 in-range=%d, CWD in-range=%d (of %d specimens)\n",
              source, sum(ok_fog), sum(ok_bio), sum(ok_cwd), n))
  for (gy in sort(unique(samples$GermYear[ok_any]))) {
    idx_fog  <- which(samples$GermYear == gy & ok_fog)
    idx_bio  <- which(samples$GermYear == gy & ok_bio)
    idx_cwd  <- which(samples$GermYear == gy & ok_cwd)
    idx_all  <- union(union(idx_fog, idx_bio), idx_cwd)
    lats  <- samples$Lat[idx_all]; longs <- samples$Long[idx_all]
    win_years <- (gy - ROLL_WINDOW + 1):gy
    cat(sprintf("    GermYear %d: window %d-%d (fog:%d, bio:%d, cwd:%d specimen(s))\n",
                gy, min(win_years), max(win_years), length(idx_fog), length(idx_bio), length(idx_cwd)))
    fog_mat <- bio12_mat <- bio15_mat <- bio18_mat <- cwd_mat <-
      matrix(NA_real_, nrow = length(idx_all), ncol = ROLL_WINDOW)
    for (i in seq_along(win_years)) {
      yr <- win_years[i]
      if (source == "modern") {
        if (yr >= fog_year_min && yr <= fog_year_max) fog_mat[, i] <- extract_fog_seasonal_13(yr, lats, longs)
      } else {
        fog_mat[, i] <- extract_kling_seasonal_15(yr, lats, longs)
      }
      pv <- extract_prism_year_13(yr, lats, longs)
      bio12_mat[, i] <- pv$BIO12; bio15_mat[, i] <- pv$BIO15; bio18_mat[, i] <- pv$BIO18
      cwd_mat[, i] <- extract_cwd_annual_15(yr, lats, longs)
    }
    fog_means   <- safe_rowMeans(fog_mat)
    bio12_means <- safe_rowMeans(bio12_mat)
    bio15_means <- safe_rowMeans(bio15_mat)
    bio18_means <- safe_rowMeans(bio18_mat)
    cwd_means   <- safe_rowMeans(cwd_mat)
    val_fog[idx_fog]     <- fog_means[match(idx_fog, idx_all)]
    val_bio12[idx_bio]   <- bio12_means[match(idx_bio, idx_all)]
    val_bio15[idx_bio]   <- bio15_means[match(idx_bio, idx_all)]
    val_bio18[idx_bio]   <- bio18_means[match(idx_bio, idx_all)]
    val_cwd[idx_cwd]     <- cwd_means[match(idx_cwd, idx_all)]
  }
  list(fog = val_fog, bio12 = val_bio12, bio15 = val_bio15, bio18 = val_bio18, cwd = val_cwd)
}

modern_fog_min   <- fog_year_min + ROLL_WINDOW - 1   # 2009
modern_fog_max   <- fog_year_max                      # 2025
modern_bio_min   <- 1927 + ROLL_WINDOW - 1
modern_bio_max   <- prism_wy_last                     # 2024
historic_fog_min <- min(kling_years_15) + ROLL_WINDOW - 1  # 1969 (20yr)
historic_fog_max <- max(kling_years_15)                     # 2014
historic_bio_min <- 1927 + ROLL_WINDOW - 1
historic_bio_max <- prism_wy_last                           # 2024 (much wider than Kling)
# CWD cache range -> same gate used for both modern and historic calls
# (unlike fog, CWD has one shared cache, not a separate modern/historic
# source), matching how BIO12/15/18's gate is already shared between them.
cwd_gate_min <- if (length(cwd_cache_years) > 0) min(cwd_cache_years) + ROLL_WINDOW - 1 else Inf
cwd_gate_max <- if (length(cwd_cache_years) > 0) max(cwd_cache_years) else -Inf

cat("\n--- modern (fog/BIO12/15/18) ---\n")
modern_vals <- compute_rolling_15(samples15, "modern", modern_fog_min, modern_fog_max,
                                   modern_bio_min, modern_bio_max, cwd_gate_min, cwd_gate_max)
cat("\n--- historic (Kling/BIO12/15/18/CWD) ---\n")
historic_vals <- compute_rolling_15(samples15, "historic", historic_fog_min, historic_fog_max,
                                     historic_bio_min, historic_bio_max, cwd_gate_min, cwd_gate_max)

samples15$fog_modern   <- modern_vals$fog
samples15$BIO12_modern <- modern_vals$bio12
samples15$BIO15_modern <- modern_vals$bio15
samples15$BIO18_modern <- modern_vals$bio18
# Not adding a "CWD_modern" column here: the modern pipeline already has its
# own established CWD_25yr (fixed WY2000-2025 mean, Section 12) -- this
# rolling-window CWD value is computed above only because compute_rolling_15
# is shared between the modern/historic calls, but only the historic side
# (which had no CWD variable at all before) actually needs a new column.
samples15$Kling_historic <- historic_vals$fog
samples15$BIO12_historic <- historic_vals$bio12
samples15$BIO15_historic <- historic_vals$bio15
samples15$BIO18_historic <- historic_vals$bio18
samples15$CWD_historic   <- historic_vals$cwd

cat(sprintf("\nNon-NA counts: fog_modern=%d BIO18_modern=%d BIO15_modern=%d | Kling_historic=%d BIO18_historic=%d BIO15_historic=%d CWD_historic=%d\n",
            sum(!is.na(samples15$fog_modern)), sum(!is.na(samples15$BIO18_modern)), sum(!is.na(samples15$BIO15_modern)),
            sum(!is.na(samples15$Kling_historic)), sum(!is.na(samples15$BIO18_historic)), sum(!is.na(samples15$BIO15_historic)),
            sum(!is.na(samples15$CWD_historic))))

# --- merge into emma_gee_results_final.csv (update existing rows by
#     SampleID, add any genuinely new ones) ---

existing15 <- read.csv(output_file, stringsAsFactors = FALSE)
existing15$SampleID <- trimws(as.character(existing15$SampleID))
existing15$Fog_historic_spliced <- NULL  # drop stale column from the now-reverted splice attempt
new_cols <- c("fog_modern", "BIO12_modern", "BIO15_modern", "BIO18_modern",
              "Kling_historic", "BIO12_historic", "BIO15_historic", "BIO18_historic",
              "CWD_historic")
for (col in new_cols) existing15[[col]] <- NA_real_

match_idx15 <- match(samples15$SampleID, existing15$SampleID)
have_match  <- !is.na(match_idx15)
for (col in new_cols) existing15[[col]][match_idx15[have_match]] <- samples15[[col]][have_match]

brand_new15 <- samples15[!have_match, ]
if (nrow(brand_new15) > 0) {
  new_rows15 <- data.frame(SampleID = brand_new15$SampleID, Lat = brand_new15$Lat, Long = brand_new15$Long,
                           Year = brand_new15$GermYear)
  for (col in setdiff(names(existing15), names(new_rows15))) new_rows15[[col]] <- NA
  new_rows15 <- new_rows15[, names(existing15)]
  for (col in new_cols) new_rows15[[col]] <- brand_new15[[col]]
  existing15 <- rbind(existing15, new_rows15)
  cat(sprintf("\n%d brand-new SampleID(s) added\n", nrow(brand_new15)))
}

write.csv(existing15, output_file, row.names = FALSE)
cat("\nRolling _modern/_historic variables written to", output_file, "\n")
