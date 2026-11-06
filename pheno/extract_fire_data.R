# extract_fire_data.R
#
# Extraction-only adaptation of the undergrad's fire.R (which mixed
# extraction, GLMs, and map-making in one file, and relied on
# hardcoded Google Drive paths plus objects — `merged`, `camap` —
# that were never defined anywhere in that file). This script keeps
# only the extraction logic (spatial join of specimen coordinates
# against fire perimeters + prescribed burns), applied to the same
# specimen set/coordinate source used throughout this project's
# other scripts (gtinatandvars.csv + inatrays.csv, full range, no
# year filter — matching historic_data.R's own scope), and writes
# the result to fire_data.txt for phenotype_data.R to read in.
#
# For each specimen, computes:
#   years_since_fire : (collection year) - (most recent prior fire's year),
#                       always >= 1 when present (positive = years elapsed
#                       since that fire); NA if never burned
#                       before collection (raw extraction keeps this as
#                       NA, not fire.R's -250 sentinel — that's a modeling
#                       choice applied downstream, not a property of the
#                       data itself)
#   burned            : "YES"/"NO"
#   fire_acreage      : GIS_ACRES of that most recent prior fire
#   fire_cause        : CAUSE of that most recent prior fire (perimeter
#                        fires only — prescribed burns get cause "0", as
#                        the undergrad's script also set)
#   n_fires_before    : count of fires (perimeter + prescribed) strictly
#                        before the collection year
#
# Output: /private/groups/meyerlab/esteiger/pheno/fire_data.txt
#         (comma-separated despite the .txt extension, matching this
#         project's other data files)

suppressMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
})

pheno_file        <- "/private/groups/meyerlab/esteiger/pheno/gtinatandvars.csv"
inatrays_file     <- "/private/groups/meyerlab/esteiger/pheno/inatrays.csv"
fire_perim_file   <- "/private/groups/meyerlab/esteiger/env/fire/California_Fire_Perimeters/California_Fire_Perimeters_(all).shp"
prescribed_file   <- "/private/groups/meyerlab/esteiger/env/fire/ds397.gdb"
output_file       <- "/private/groups/meyerlab/esteiger/pheno/fire_data.txt"

normalize_locality <- function(x) {
  # inatrays.csv has some rows with "VAFB_Lion_head" (lowercase h) instead
  # of "VAFB_Lion_Head" — same inconsistency plot_env_time.R already
  # works around; fixed here before stripping underscores so it doesn't
  # silently split into two localities.
  x <- gsub("_head$", "_Head", x)
  x <- gsub("_", "", x)
  as.character(ifelse(is.na(x) | x == "", NA_character_, x))
}

# ============================================================
# 1. Assemble the same specimen set used elsewhere (full range, no year
#    filter — matches historic_data.R's scope)
# ============================================================

pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
pheno$recordNumber_CollectorNumber <- trimws(as.character(pheno$recordNumber_CollectorNumber))
pheno$year <- suppressWarnings(as.numeric(pheno$year))

samples <- pheno %>%
  filter(!is.na(year), !is.na(Lat._dd.dddd), !is.na(Lon_dd.dddd)) %>%
  transmute(SampleID = recordNumber_CollectorNumber, year,
            Locality = normalize_locality(Locality),
            Lat = Lat._dd.dddd, Long = Lon_dd.dddd) %>%
  distinct()

inat <- read.csv(inatrays_file, stringsAsFactors = FALSE)
inat$inat_Specimen <- trimws(as.character(inat$inat_Specimen))
inat_ok <- inat %>%
  filter(!is.na(Year), !is.na(latitude), !is.na(longitude),
         !(inat_Specimen %in% samples$SampleID)) %>%
  transmute(SampleID = inat_Specimen, year = Year,
            Locality = normalize_locality(Locality),
            Lat = latitude, Long = longitude) %>%
  distinct()

samples <- bind_rows(samples, inat_ok) %>% distinct(SampleID, .keep_all = TRUE)
cat(sprintf("Specimens with valid coordinates + year: %d\n", nrow(samples)))

proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
sf_points <- st_as_sf(samples, coords = c("Long", "Lat"), crs = proj, remove = FALSE)

# ============================================================
# 2. Load fire perimeters + prescribed burns, restricted to the study
#    area bounding box (same box fire.R used)
# ============================================================

bbox <- c(xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9)

cat("Reading fire perimeters shapefile...\n")
fire_perim <- st_read(fire_perim_file, quiet = TRUE) %>%
  st_transform(crs = proj) %>%
  st_make_valid() %>%
  st_crop(bbox)

cat("Reading prescribed-burns geodatabase (slow — large multi-polygon file)...\n")
prescribed <- st_read(prescribed_file, quiet = TRUE) %>%
  st_transform(crs = proj) %>%
  st_cast("MULTIPOLYGON") %>%
  st_make_valid() %>%
  st_crop(bbox)

cat(sprintf("Fire perimeters in study area: %d | Prescribed burns in study area: %d\n",
            nrow(fire_perim), nrow(prescribed)))

# ============================================================
# 3. Spatial join: which specimens fall inside which fire polygons
# ============================================================

perim_hits <- st_join(fire_perim, sf_points) %>%
  st_drop_geometry() %>%
  filter(!is.na(SampleID)) %>%
  transmute(SampleID, fire_year = as.numeric(YEAR_), fire_date = as.character(ALARM_DATE),
            acres = as.numeric(GIS_ACRES), cause = as.character(CAUSE))

presc_hits <- st_join(prescribed, sf_points) %>%
  st_drop_geometry() %>%
  filter(!is.na(SampleID)) %>%
  transmute(SampleID, fire_year = as.numeric(YEAR_), fire_date = as.character(START_DATE),
            acres = as.numeric(GIS_ACRES), cause = "0")

all_hits <- bind_rows(perim_hits, presc_hits)
cat(sprintf("Specimen-fire intersections found: %d (%d unique specimens)\n",
            nrow(all_hits), length(unique(all_hits$SampleID))))

# ============================================================
# 4. Per-specimen summary: most recent fire strictly before the
#    specimen's own collection year
# ============================================================

years_since_fire <- rep(NA_real_, nrow(samples))
fire_acreage     <- rep(NA_real_, nrow(samples))
fire_cause       <- rep(NA_character_, nrow(samples))
n_fires_before   <- rep(0L, nrow(samples))

for (i in seq_len(nrow(samples))) {
  sid <- samples$SampleID[i]
  cyear <- samples$year[i]
  h <- all_hits %>% filter(SampleID == sid, !is.na(fire_year))
  prior <- h %>% filter(fire_year < cyear)
  n_fires_before[i] <- nrow(prior)
  if (nrow(prior) > 0) {
    best <- prior %>% arrange(desc(fire_year)) %>% slice(1)
    years_since_fire[i] <- cyear - best$fire_year
    fire_acreage[i]     <- best$acres
    fire_cause[i]        <- best$cause
  }
}

fire_data <- samples %>%
  mutate(
    years_since_fire = years_since_fire,
    burned           = ifelse(!is.na(years_since_fire), "YES", "NO"),
    fire_acreage     = fire_acreage,
    fire_cause       = fire_cause,
    n_fires_before   = n_fires_before
  ) %>%
  select(SampleID, Locality, year, Lat, Long, years_since_fire, burned, fire_acreage, fire_cause, n_fires_before)

cat(sprintf("\nBurned before collection: %d of %d specimens\n",
            sum(fire_data$burned == "YES"), nrow(fire_data)))
cat("By locality:\n")
print(fire_data %>% count(Locality, burned))

write.csv(fire_data, output_file, row.names = FALSE)
cat(sprintf("\nFire data written to %s\n", output_file))
