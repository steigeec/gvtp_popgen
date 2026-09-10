#!/usr/bin/env Rscript
# =============================================================================
# explore_env_var.R
# =============================================================================

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")
setwd("/private/groups/meyerlab/esteiger")
outdir <- "/private/groups/meyerlab/esteiger/selection"

# ==============================================================================
# 1. READ AND MERGE
# ==============================================================================
env   <- read.table("env_var.txt",    header = TRUE, sep = "\t")
field <- read.table("field_data.txt", header = TRUE, sep = "\t")
edna  <- read.csv("eDNA_env_data.csv", header = TRUE)

env$SampleID   <- as.character(env$SampleID)
field$SampleID <- as.character(field$SampleID)
edna$SampleID  <- as.character(edna$SampleID)

merged <- merge(env,
                field[, c("SampleID", "Clusters", "index", "OHANA_order")],
                by = "SampleID", all.x = TRUE)
merged <- merge(merged, edna, by = "SampleID", all.x = TRUE)

cat(sprintf("Merged dataset: %d samples, %d columns\n\n", nrow(merged), ncol(merged)))

# aspect is circular (degrees, 0-360 wraps around) so Pearson/Spearman on the
# raw value is not meaningful. Decompose into northness/eastness (unit-circle
# projection) — two ordinary linear variables that losslessly encode the same
# direction information and can be used in correlation/PCA/GEMMA directly.
merged$aspect_northness <- cos(merged$aspect * pi / 180)
merged$aspect_eastness  <- sin(merged$aspect * pi / 180)

# ==============================================================================
# 2. CONFIRMED YES VARIABLES (explicit user decisions)
# ==============================================================================
confirmed_yes <- c(
  # Fog / precip (sampling year)
  "fog", "BCM_30y",
  # WorldClim long-term climate
  "BIO12", "BIO15", "BIO18",
  # Terrain
  "elevation", "slope", "aspect_northness", "aspect_eastness",
  # Soil (10 cm only)
  "sand_10cm", "PercVolSoilH20_10cm", "bulkDensity_10cm",
  # Land cover
  "fapar_annual_fapar"
)

# ==============================================================================
# 3. CONFIRMED NO — variables explicitly excluded from analysis
# ==============================================================================
confirmed_no <- c(
  # the germination ear variables"
  "fog_germ_yr", "BIO12_germ_yr", "BIO15_germ_yr", "BIO18_germ_yr", "BCM_germ_yr",
  # CHILI — excluded (TODO: confirm rationale)
  "CHILI",
  # textureClass_10cm — dataset only has 2 valid texture classes (4 and 7;
  # -999999 is a missing-data sentinel, not a real class) — too little
  # variation to be a useful covariate
  "textureClass_10cm",
  # dynamic_world_change_dataset_land_cover_fractional_change — land-cover
  # change / habitat loss is not expected to be a direct driver of selection
  # in this system, so excluded even though it is a valid numeric variable
  "dynamic_world_change_dataset_land_cover_fractional_change",
  # TerraClimate (~4.6 km resolution)
  "terraclimate_climate_water_deficit_mm_30y",
  "terraclimate_soil_moisture_mm_30y",
  # Copernicus fractional cover — not ecologically informative for this project
  "copernicus_fractional_cover_tree_cover_fraction",
  "copernicus_fractional_cover_shrub_cover_fraction",
  "terraclimate_snow_water_equivalent_mm_30y",
  "copernicus_fractional_cover_grass_cover_fraction",
  "copernicus_fractional_cover_bare_cover_fraction",
  "copernicus_fractional_cover_crops_cover_fraction",
  "copernicus_fractional_cover_urban_cover_fraction",
  "copernicus_fractional_cover_water_permanent_cover_fraction",
  "copernicus_fractional_cover_water_seasonal_cover_fraction",
  "copernicus_fractional_cover_moss_lichen_cover_fraction",
  # Hansen forest cover/loss — not of interest
  "hansen_forest_forest_cover_2000_percent",
  "hansen_forest_forest_loss_binary_2001_2024",
  "hansen_forest_forest_loss_rate_pct",
  # JRC surface water — not of interest
  "jrc_surface_water_max_water_extent_binary", 
  # Any 10-year variables (we're interested in the 30-year timespan)
  "era5_land_mean_net_solar_radiation_w_m2_10y",
  "terraclimate_downward_surface_shortwave_radiation_w_m2_10y",
  "terraclimate_actual_evapotranspiration_mm_10y",
  "terraclimate_palmer_drought_severity_index_10y",
  "terraclimate_potential_evapotranspiration_mm_10y",
  "terraclimate_vapor_pressure_deficit_kpa_10y",
  "terraclimate_vapor_pressure_kpa_10y",
  "terraclimate_wind_speed_m_s_10y",
  "terraclimate_snow_water_equivalent_mm_10y",
  "terraclimate_aridity_index_10y",
  "terraclimate_soil_moisture_seasonality_cv_10y",
  "terraclimate_precip_annual_cv_10y",
  "terraclimate_severe_drought_month_count_10y",
  "terraclimate_driest_year_precip_mm_10y",
  # Additional TerraClimate 10y variables (same rationale as above)
  "terraclimate_annual_mean_temp_10y",
  "terraclimate_annual_precip_10y",
  "terraclimate_isothermality_10y",
  "terraclimate_max_temp_warmest_month_10y",
  "terraclimate_mean_diurnal_range_10y",
  "terraclimate_mean_temp_coldest_quarter_10y",
  "terraclimate_mean_temp_driest_quarter_10y",
  "terraclimate_mean_temp_warmest_quarter_10y",
  "terraclimate_mean_temp_wettest_quarter_10y",
  "terraclimate_min_temp_coldest_month_10y",
  "terraclimate_precip_coldest_quarter_10y",
  "terraclimate_precip_driest_month_10y",
  "terraclimate_precip_driest_quarter_10y",
  "terraclimate_precip_seasonality_10y",
  "terraclimate_precip_warmest_quarter_10y",
  "terraclimate_precip_wettest_month_10y",
  "terraclimate_precip_wettest_quarter_10y",
  "terraclimate_temp_annual_range_10y",
  "terraclimate_temp_seasonality_10y",
  "terraclimate_climate_water_deficit_mm_10y",
  "terraclimate_max_temp_c_10y",
  "terraclimate_min_temp_c_10y",
  "terraclimate_precipitation_accumulation_mm_10y",
  "terraclimate_runoff_mm_10y",
  "terraclimate_soil_moisture_mm_10y",
  "terraclimate_gdd_10c_10y",
  "terraclimate_gdd_5c_10y",
  "terraclimate_growing_season_length_10c_10y",
  "terraclimate_growing_season_length_5c_10y",
  "terraclimate_hottest_year_mean_temp_c_10y",
  # Watershed — not ecologically relevant for this project
  "watershed_watershed_id",
  "watershed_watershed_pfaf_id",
  "watershed_watershed_upstream_area_km2",
  "watershed_watershed_sub_area_km2",
  "watershed_watershed_dist_to_main_outlet_km",
  "watershed_watershed_endorheic_indicator",
  "watershed_watershed_coastal_indicator",
  "watershed_watershed_river_order",
  # CHIRPS precipitation metrics
  "chirps_heavy_precip_days_count_30y",
  "chirps_rainfall_intensity_sdii_mm_day_30y",
  "chirps_heavy_precip_days_count_10y",
  "chirps_rainfall_intensity_sdii_mm_day_10y",
  # MODIS land-surface temperature
  "modis_annual_max_day_lst_c_3y",
  "modis_annual_mean_day_lst_c_3y",
  "modis_annual_max_day_lst_c_1y",
  "modis_annual_mean_day_lst_c_1y",
  # MODIS vegetation indices
  "modis_vegetation_indices_evi_seasonality_cv_1y",
  "modis_vegetation_indices_growing_season_days_1y",
  "modis_vegetation_indices_ndvi_annual_amplitude_1y",
  "modis_vegetation_indices_peak_greenness_doy_1y",
  "modis_vegetation_indices_senescence_doy_1y",
  "modis_vegetation_indices_evi_seasonality_cv_3y",
  "modis_vegetation_indices_growing_season_days_3y",
  "modis_vegetation_indices_ndvi_annual_amplitude_3y",
  "modis_vegetation_indices_peak_greenness_doy_3y",
  "modis_vegetation_indices_senescence_doy_3y",
  "modis_vegetation_indices_evi_seasonality_cv_10y",
  "modis_vegetation_indices_growing_season_days_10y",
  "modis_vegetation_indices_ndvi_annual_amplitude_10y",
  "modis_vegetation_indices_peak_greenness_doy_10y",
  "modis_vegetation_indices_senescence_doy_10y",
  # Sentinel-5P atmospheric chemistry
  "sentinel5p_mean_annual_no2_mol_m2_1y",
  "sentinel5p_mean_annual_no2_mol_m2_3y",
  "sentinel5p_mean_annual_co_mol_m2_1y",
  "sentinel5p_mean_annual_co_mol_m2_3y",
  "sentinel5p_mean_annual_so2_mol_m2_1y",
  "sentinel5p_mean_annual_so2_mol_m2_3y",
  "sentinel5p_mean_annual_o3_mol_m2_1y",
  "sentinel5p_mean_annual_o3_mol_m2_3y",
  "sentinel5p_mean_annual_hcho_mol_m2_1y",
  "sentinel5p_mean_annual_hcho_mol_m2_3y",
  "sentinel5p_mean_annual_aerosol_index_1y",
  "sentinel5p_mean_annual_aerosol_index_3y",
  # Nighttime lights
  "nighttime_lights_max_annual_nightlights_nw_cm2_sr_1y",
  "nighttime_lights_mean_annual_nightlights_nw_cm2_sr_1y",
  "nighttime_lights_mean_ephemeral_nightlights_1y",
  "nighttime_lights_mean_persistent_nightlights_1y",
  "nighttime_lights_max_annual_nightlights_nw_cm2_sr_3y",
  "nighttime_lights_mean_annual_nightlights_nw_cm2_sr_3y",
  "nighttime_lights_mean_ephemeral_nightlights_3y",
  "nighttime_lights_mean_persistent_nightlights_3y",
  # Metadata / identifiers — not analytical variables
  "protected_areas_protected_area_name",
  "protected_areas_protected_area_iucn_cat",
  "resolve_ecoregions_ecoregion_name_2017",
  "resolve_ecoregions_biogeographic_realm_2017",
  # Terrain duplicates — elevation, slope, aspect already in confirmed_yes
  "terrain_elevation_m",
  "terrain_slope_deg",
  "terrain_aspect_deg",
  "terrain_topo_heterogeneity_stddev",
  # Soil duplicates — corresponding variables already in confirmed_yes
  "soil_soil_sand_percent_0_30cm", "soil_soil_sand_percent_0cm",
  "soil_soil_bulk_density_cg_cm3_0_30cm", "soil_soil_bulk_density_cg_cm3_0cm",
  "soil_soil_water_content_v_percent_0_30cm", "soil_soil_water_content_v_percent_0cm",
  "soil_soil_texture_class_0cm",
  "soil_soil_great_group",
  # New soil properties — 0 cm depth excluded (0–30 cm retained as candidate)
  "soil_soil_ph_h2o_0cm",
  "soil_soil_organic_carbon_g_kg_0cm",
  "soil_soil_clay_percent_0cm",
  # New Copernicus columns — not ecologically informative for this project
  "copernicus_fractional_cover_snow_ice_cover_fraction",
  "copernicus_fractional_cover_class_probability",
  "copernicus_fractional_cover_land_cover_class",
  "copernicus_fractional_cover_forest_type",
  # Land cover / biome categorical — not informative as PCA/RF axes
  "esa_worldcover_dataset_land_cover_class_2020",
  "potential_biomes_potential_biome",
  # TerraClimate temperature (30y) — not a primary driver for this project
  "terraclimate_annual_mean_temp_30y",
  "terraclimate_max_temp_c_30y",
  "terraclimate_min_temp_c_30y",
  "terraclimate_max_temp_warmest_month_30y",
  "terraclimate_min_temp_coldest_month_30y",
  "terraclimate_mean_diurnal_range_30y",
  "terraclimate_isothermality_30y",
  "terraclimate_temp_annual_range_30y",
  "terraclimate_temp_seasonality_30y",
  "terraclimate_mean_temp_coldest_quarter_30y",
  "terraclimate_mean_temp_driest_quarter_30y",
  "terraclimate_mean_temp_warmest_quarter_30y",
  "terraclimate_mean_temp_wettest_quarter_30y",
  "terraclimate_hottest_year_mean_temp_c_30y",
  "landsat_peak_season_surface_temp_c",
  # TerraClimate precipitation (30y) — covered by BIO12/BIO15 and existing candidates
  "terraclimate_annual_precip_30y",
  "terraclimate_precipitation_accumulation_mm_30y",
  "terraclimate_precip_seasonality_30y",
  "terraclimate_precip_coldest_quarter_30y",
  "terraclimate_precip_driest_month_30y",
  "terraclimate_precip_driest_quarter_30y",
  "terraclimate_precip_warmest_quarter_30y",
  "terraclimate_precip_wettest_month_30y",
  "terraclimate_precip_wettest_quarter_30y",
  "terraclimate_runoff_mm_30y",
  # TerraClimate growing season / phenology (30y) — not primary interest
  "terraclimate_gdd_10c_30y",
  "terraclimate_gdd_5c_30y",
  "terraclimate_growing_season_length_10c_30y",
  "terraclimate_growing_season_length_5c_30y"
)

# ==============================================================================
# 4. CANDIDATE POOL — variables still to be decided
#    (excludes confirmed YES and confirmed NO)
# ==============================================================================
candidates <- c(
  # Radiation / energy
  "era5_land_mean_net_solar_radiation_w_m2_30y",
  "terraclimate_downward_surface_shortwave_radiation_w_m2_30y",
  # TerraClimate climate (undecided)
  "terraclimate_actual_evapotranspiration_mm_30y",
  "terraclimate_palmer_drought_severity_index_30y",
  "terraclimate_potential_evapotranspiration_mm_30y",
  "terraclimate_vapor_pressure_deficit_kpa_30y",
  "terraclimate_vapor_pressure_kpa_30y",
  "terraclimate_wind_speed_m_s_30y",
  "terraclimate_aridity_index_30y",
  "terraclimate_soil_moisture_seasonality_cv_30y",
  "terraclimate_precip_annual_cv_30y",
  "terraclimate_severe_drought_month_count_30y",
  "terraclimate_driest_year_precip_mm_30y",
  # New SoilGrids soil properties (250 m, 0–30 cm depth)
  "soil_soil_ph_h2o_0_30cm",
  "soil_soil_organic_carbon_g_kg_0_30cm",
  "soil_soil_clay_percent_0_30cm"
)
candidates <- candidates[candidates %in% names(merged)]

# ==============================================================================
# 5. STRAUSS SUBSETTING
#    Strauss is heavily over-represented; subsample to 10 random individuals
#    before missingness filtering so that missingness reflects the analysis
#    dataset, not the Strauss-dominated full dataset.
#    (Full Strauss data retained in `merged` for the final export.)
# ==============================================================================
library(ggplot2)
library(randomForest)

set.seed(42)
cluster_counts <- table(merged$Clusters[!is.na(merged$Clusters)])
n_strauss_sub  <- max(cluster_counts[names(cluster_counts) != "Strauss"])
cat(sprintf("Largest non-Strauss cluster: n = %d — subsampling Strauss to match\n",
            n_strauss_sub))
strauss_ids  <- merged$SampleID[!is.na(merged$Clusters) & merged$Clusters == "Strauss"]
strauss_keep <- sample(strauss_ids, n_strauss_sub)
analysis_data <- merged[!is.na(merged$Clusters) &
                        (merged$Clusters != "Strauss" |
                         merged$SampleID %in% strauss_keep), ]
n_strauss_kept <- sum(analysis_data$Clusters == "Strauss", na.rm = TRUE)
cat(sprintf("Strauss subsampled: %d individuals kept (of %d total)\n",
            n_strauss_kept, length(strauss_ids)))
cat(sprintf("Analysis dataset: %d samples\n\n", nrow(analysis_data)))

# ==============================================================================
# 6. FILTER 1 — missingness (> 10% missing → drop)
#    Applied to analysis_data; checks both confirmed_yes and candidate variables.
# ==============================================================================
n <- nrow(analysis_data)

# Check confirmed_yes variables for missingness in analysis subset
all_check <- c(confirmed_yes, candidates)
all_check  <- all_check[all_check %in% names(analysis_data)]
pct_miss   <- sapply(all_check, function(v) 100 * sum(is.na(analysis_data[[v]])) / n)

# Confirmed YES variables — report missingness for all, flag those > threshold
cat("Confirmed-YES variable missingness in analysis subset:\n")
for (v in confirmed_yes[confirmed_yes %in% names(analysis_data)]) {
  flag <- if (pct_miss[v] > 10) "  *** EXCEEDS 10% — add to dropped_for_missingness" else ""
  cat(sprintf("  %-55s %.1f%%%s\n", v, pct_miss[v], flag))
}

# Candidates that fail missingness
fail_miss <- names(pct_miss[pct_miss > 10 & names(pct_miss) %in% candidates])
pass_miss <- candidates[!candidates %in% fail_miss]

cat(sprintf("\nCandidate missingness filter (>10%% missing → drop):\n"))
cat(sprintf("  Dropped: %d  |  Remaining: %d\n", length(fail_miss), length(pass_miss)))
if (length(fail_miss)) cat("  Dropped:", paste(fail_miss, collapse=", "), "\n")

# ==============================================================================
# 7. FILTER 2 — spatial resolution (> 1 km → drop)
#    Native resolutions (metres) per data source:
#      SRTM / terrain derived : 30 m
#      SoilGrids              : 250 m
#      Copernicus LC          : 100 m
#      Hansen Global Forest   : 30 m
#      JRC Surface Water      : 30 m
#      Dynamic World          : 10 m
#      MODIS (FAPAR, LST)     : 500–1000 m
#      WorldClim BIO          : ~1000 m (30 arc-sec)
#      MODIS FLCC (fog)       : 1000 m
#      TerraClimate           : ~4638 m (1/24°)  ← EXCEEDS 1 km
#      ERA5-Land              : ~11132 m (0.1°)  ← EXCEEDS 1 km
# ==============================================================================
resolution_m <- c(
  era5_land_mean_net_solar_radiation_w_m2_30y              = 11132,
  era5_land_mean_net_solar_radiation_w_m2_10y              = 11132,
  terraclimate_downward_surface_shortwave_radiation_w_m2_30y = 4638,
  terraclimate_downward_surface_shortwave_radiation_w_m2_10y = 4638,
  terraclimate_actual_evapotranspiration_mm_30y            = 4638,
  terraclimate_actual_evapotranspiration_mm_10y            = 4638,
  terraclimate_palmer_drought_severity_index_30y           = 4638,
  terraclimate_palmer_drought_severity_index_10y           = 4638,
  terraclimate_potential_evapotranspiration_mm_30y         = 4638,
  terraclimate_potential_evapotranspiration_mm_10y         = 4638,
  terraclimate_vapor_pressure_deficit_kpa_30y              = 4638,
  terraclimate_vapor_pressure_deficit_kpa_10y              = 4638,
  terraclimate_vapor_pressure_kpa_30y                      = 4638,
  terraclimate_vapor_pressure_kpa_10y                      = 4638,
  terraclimate_wind_speed_m_s_30y                          = 4638,
  terraclimate_wind_speed_m_s_10y                          = 4638,
  terraclimate_snow_water_equivalent_mm_30y                = 4638,
  terraclimate_snow_water_equivalent_mm_10y                = 4638,
  terraclimate_aridity_index_30y                           = 4638,
  terraclimate_aridity_index_10y                           = 4638,
  terraclimate_soil_moisture_seasonality_cv_30y            = 4638,
  terraclimate_soil_moisture_seasonality_cv_10y            = 4638,
  terraclimate_precip_annual_cv_30y                        = 4638,
  terraclimate_precip_annual_cv_10y                        = 4638,
  terraclimate_severe_drought_month_count_30y              = 4638,
  terraclimate_severe_drought_month_count_10y              = 4638,
  terraclimate_driest_year_precip_mm_30y                   = 4638,
  terraclimate_driest_year_precip_mm_10y                   = 4638,
  # SoilGrids — new soil property candidates
  soil_soil_ph_h2o_0_30cm                                  = 250,
  soil_soil_organic_carbon_g_kg_0_30cm                     = 250,
  soil_soil_clay_percent_0_30cm                            = 250
)

res_pass_miss <- resolution_m[pass_miss]
fail_res <- names(res_pass_miss)[res_pass_miss > 1000 & !is.na(res_pass_miss)]
pass_both <- pass_miss[!pass_miss %in% fail_res]

cat(sprintf("\nResolution filter (>1 km → drop):\n"))
cat(sprintf("  Dropped: %d  |  Remaining: %d\n", length(fail_res), length(pass_both)))
if (length(fail_res)) cat("  Dropped:", paste(fail_res, collapse="\n           "), "\n")

# ==============================================================================
# 8. SURVIVING CANDIDATES (after both filters)
# ==============================================================================
cat("\n=== Surviving candidates after missingness + resolution filters ===\n")
surv <- data.frame(
  variable   = pass_both,
  res_m      = resolution_m[pass_both],
  pct_miss   = round(pct_miss[pass_both], 1),
  row.names  = NULL
)
print(surv, row.names = FALSE)

COLLINEAR_THRESH <- 0.7

# No categorical/circular variables remain — aspect was replaced with the
# numeric aspect_northness/aspect_eastness decomposition (see section 1).
excl_nonnumeric <- character(0)

# ==============================================================================
# 9. PRE-WEEDING COLLINEARITY HEATMAP
#    All confirmed_yes + surviving candidates, before any collinear removal.
#    Review the flagged pairs and heatmap to populate remove_collinear in
#    section 10.
# ==============================================================================
prefilter_pool <- c(confirmed_yes, pass_both)
prefilter_num  <- prefilter_pool[!prefilter_pool %in% excl_nonnumeric]
prefilter_num  <- prefilter_num[prefilter_num %in% names(analysis_data)]

cor_mat_pre <- cor(analysis_data[, prefilter_num], method = "spearman",
                   use = "pairwise.complete.obs")

pairs_idx_pre <- which(upper.tri(cor_mat_pre), arr.ind = TRUE)
cor_pairs_pre <- data.frame(
  var1 = rownames(cor_mat_pre)[pairs_idx_pre[, 1]],
  var2 = colnames(cor_mat_pre)[pairs_idx_pre[, 2]],
  rho  = round(cor_mat_pre[pairs_idx_pre], 3)
)
cor_pairs_pre <- cor_pairs_pre[order(-abs(cor_pairs_pre$rho)), ]
flagged_pre   <- cor_pairs_pre[abs(cor_pairs_pre$rho) >= COLLINEAR_THRESH, ]

cat(sprintf("\n=== Pre-weeding collinearity (Spearman |rho| >= %.1f) ===\n", COLLINEAR_THRESH))
cat(sprintf("Flagged pairs: %d\n\n", nrow(flagged_pre)))
print(flagged_pre, row.names = FALSE)

hc_ord_pre     <- hclust(as.dist(1 - cor_mat_pre))$order
var_levels_pre <- rownames(cor_mat_pre)[hc_ord_pre]

cor_df_pre <- data.frame(
  var1 = rep(rownames(cor_mat_pre), ncol(cor_mat_pre)),
  var2 = rep(colnames(cor_mat_pre), each = nrow(cor_mat_pre)),
  rho  = as.vector(cor_mat_pre)
)
cor_df_pre$var1 <- factor(cor_df_pre$var1, levels = var_levels_pre)
cor_df_pre$var2 <- factor(cor_df_pre$var2, levels = var_levels_pre)

jpeg(file.path(outdir, "env_correlation_heatmap_preweeding.jpg"), width = 1600, height = 1400)
print(
  ggplot(cor_df_pre, aes(x = var2, y = var1, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title   = element_blank()) +
    ggtitle("Pre-weeding Spearman correlations — confirmed_yes + candidates (Strauss n=10)")
)
dev.off()
cat(sprintf("Pre-weeding heatmap written to %s\n", file.path(outdir, "env_correlation_heatmap_preweeding.jpg")))

# ==============================================================================
# 10. REMOVE COLLINEAR VARIABLES
#     Populate based on flagged pairs from the pre-weeding heatmap above.
# ==============================================================================
remove_collinear <- c("bulkDensity_10cm", # remove instead of fog (-0.70)
			"terraclimate_wind_speed_m_s_30y", #remove instead of BIO15, precip seasonality (-0.71)
			"elevation", # remove instead of BIO18, precip in warmest quarter (0.70)
			"soil_soil_clay_percent_0_30cm", # remove instead of fog (-0.77)
			"soil_soil_organic_carbon_g_kg_0_30cm", #remove instead of fapar (0.83)
			"sand_10cm", # remove instead of fog (0.83)
			"slope", # remove instead of BIO12 (0.645, below the 0.7 auto-flag threshold but visibly collinear on the heatmap)
			"BIO12", # remove instead of BIO18 (0.697, below the 0.7 auto-flag threshold but visibly collinear on the heatmap)
			"PercVolSoilH20_10cm", # remove instead of fog (-0.655, below the 0.7 auto-flag threshold but visibly collinear on the heatmap)
			"soil_soil_ph_h2o_0_30cm") # remove instead of BCM_30y (0.640) and BIO18 (-0.577), both below the 0.7 auto-flag threshold but visibly collinear on the heatmap
			# fapar_annual_fapar retained despite collinearity with BIO15 (-0.541) per user decision

all_vars <- c(confirmed_yes, pass_both)
env_vars <- all_vars[!all_vars %in% remove_collinear]
cat(sprintf("\nAfter removing collinear variables: %d variables retained\n", length(env_vars)))
cat(paste(" ", env_vars, collapse = "\n"), "\n")

# ==============================================================================
# 9b. VARIABLES DROPPED FOR MISSINGNESS (> 10% in analysis subset)
#     Populated after running sections 6–8 on the subsampled dataset.
#     Currently: no confirmed-YES variables exceeded the 10% threshold.
# ==============================================================================
dropped_for_missingness <- c(
  # none
)
env_vars <- env_vars[!env_vars %in% dropped_for_missingness]

# ==============================================================================
# 10. ANALYSIS VARIABLES FOR COLLINEARITY / PCA 
#     fog_germ_yr and BIO12_germ_yr excluded; fog and BIO12 used instead.
#     The _germ_yr columns are retained in the exported file (section 13).
# ==============================================================================
# Variables for collinearity / PCA / RF: exclude _germ_yr versions
analysis_vars <- env_vars[!env_vars %in% c("fog_germ_yr", "BIO12_germ_yr")]

# Numeric subset for correlation / PCA
# No exclusions needed — aspect_northness/aspect_eastness are ordinary
# numeric variables (see the northness/eastness decomposition in section 1).
pca_exclude <- character(0)
analysis_num <- analysis_vars[!analysis_vars %in% pca_exclude]
analysis_num <- analysis_num[analysis_num %in% names(analysis_data)]

# ==============================================================================
# 11. COLLINEARITY — Spearman correlations (analysis subset, no _germ_yr)
#     Pairs with |rho| >= 0.7 flagged.
# ==============================================================================
COLLINEAR_THRESH <- 0.7

cor_mat <- cor(analysis_data[, analysis_num], method = "spearman",
               use = "pairwise.complete.obs")

pairs_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
cor_pairs <- data.frame(
  var1 = rownames(cor_mat)[pairs_idx[, 1]],
  var2 = colnames(cor_mat)[pairs_idx[, 2]],
  rho  = round(cor_mat[pairs_idx], 3)
)
cor_pairs <- cor_pairs[order(-abs(cor_pairs$rho)), ]
flagged_pairs <- cor_pairs[abs(cor_pairs$rho) >= COLLINEAR_THRESH, ]

cat(sprintf("\n=== Collinearity (Spearman |rho| >= %.1f) ===\n", COLLINEAR_THRESH))
cat(sprintf("Flagged pairs: %d\n\n", nrow(flagged_pairs)))
print(flagged_pairs, row.names = FALSE)

# Correlation heatmap
hc_ord     <- hclust(as.dist(1 - cor_mat))$order
var_levels <- rownames(cor_mat)[hc_ord]

# Display labels for the plot only (underlying variable names elsewhere in
# the script, e.g. env_vars/analysis_num, are untouched)
plot_labels <- c(
  BCM_30y            = "CWD",
  fog                = "FLCC",
  fapar_annual_fapar = "FAPAR",
  BIO18              = "BIO18",
  aspect_northness   = "northness",
  aspect_eastness    = "eastness",
  BIO15              = "BIO15"
)

cor_df <- data.frame(
  var1 = rep(rownames(cor_mat), ncol(cor_mat)),
  var2 = rep(colnames(cor_mat), each = nrow(cor_mat)),
  rho  = as.vector(cor_mat)
)
cor_df$var1 <- factor(plot_labels[as.character(cor_df$var1)],
                      levels = plot_labels[var_levels])
cor_df$var2 <- factor(plot_labels[as.character(cor_df$var2)],
                      levels = plot_labels[var_levels])

# Sized and font-scaled for print (e.g. a single figure column in a paper)
# rather than on-screen viewing: physical inches + 300 dpi so the figure is
# still legible at a small placed size, instead of a huge low-density raster
# that goes illegible when shrunk down.
jpeg(file.path(outdir, "env_correlation_heatmap.jpg"),
     width = 5, height = 4.5, units = "in", res = 300)
print(
  ggplot(cor_df, aes(x = var2, y = var1, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title   = element_blank())
)
dev.off()
cat(sprintf("Heatmap written to %s\n", file.path(outdir, "env_correlation_heatmap.jpg")))

# ==============================================================================
# 12. PCA (analysis subset, numeric analysis vars)
# ==============================================================================
pca_data <- analysis_data[, analysis_num]
pca_data <- pca_data[complete.cases(pca_data), ]
pca_sids <- analysis_data$SampleID[complete.cases(analysis_data[, analysis_num])]
cat(sprintf("\nPCA: %d variables, %d complete-case samples\n",
            ncol(pca_data), nrow(pca_data)))

pca_res <- prcomp(pca_data, center = TRUE, scale. = TRUE)

var_exp <- summary(pca_res)$importance
cat("\nVariance explained:\n")
print(round(var_exp[, 1:min(10, ncol(var_exp))], 3))

n_pc     <- min(10, ncol(pca_res$rotation))
cat("\nLoadings (PC1–PC", n_pc, "):\n", sep = "")
print(round(pca_res$rotation[, 1:n_pc], 3))

pct_var  <- 100 * pca_res$sdev^2 / sum(pca_res$sdev^2)
scree_df <- data.frame(PC = seq_along(pct_var), pct_var = pct_var,
                       cumvar = cumsum(pct_var))

jpeg(file.path(outdir, "env_pca_scree.jpg"), width = 900, height = 700)
print(
  ggplot(scree_df, aes(x = PC, y = pct_var)) +
    geom_col(fill = "steelblue") +
    geom_line(aes(y = cumvar)) +
    geom_point(aes(y = cumvar)) +
    geom_hline(yintercept = 80, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = seq_len(nrow(scree_df))) +
    labs(x = "PC", y = "Variance explained (%)") +
    theme_bw(base_size = 14)
)
dev.off()
cat(sprintf("Scree plot written to %s\n", file.path(outdir, "env_pca_scree.jpg")))

scores <- as.data.frame(pca_res$x)
scores$SampleID <- pca_sids
scores <- merge(scores, merged[, c("SampleID", "Clusters")], by = "SampleID")

arrows_df  <- as.data.frame(pca_res$rotation[, 1:2])
arrows_df$var <- plot_labels[rownames(arrows_df)]
scale_arrow <- 4

popcolors <- c("CDFW"            = "#e69f00",
               "CDFW_MR"         = "#e69f00",
               "Strauss"         = "#009e73",
               "Point_Conception" = "#0072b2",
               "VAFB_Lion_Head"  = "#f0e442",
               "VAFB_Arguello"   = "#d55e00",
               "Dangermond"      = "#56b4e9",
               "Between"         = "black")

jpeg(file.path(outdir, "env_pca_biplot.jpg"), width = 1100, height = 1000)
print(
  ggplot(scores, aes(x = PC1, y = PC2, colour = Clusters)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_segment(data = arrows_df,
                 aes(x = 0, y = 0,
                     xend = PC1 * scale_arrow, yend = PC2 * scale_arrow),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey30") +
    geom_text(data = arrows_df,
              aes(x = PC1 * scale_arrow * 1.1, y = PC2 * scale_arrow * 1.1,
                  label = var),
              inherit.aes = FALSE, size = 12, colour = "grey20") +
    scale_colour_manual(values = popcolors) +
    labs(x = sprintf("PC1 (%.1f%%)", pct_var[1]),
         y = sprintf("PC2 (%.1f%%)", pct_var[2])) +
    theme_bw(base_size = 70) +
    theme(legend.position = "none")
)
dev.off()
cat(sprintf("Biplot written to %s\n", file.path(outdir, "env_pca_biplot.jpg")))

# ==============================================================================
# 13. FINAL SET — second collinearity + PCA pass excluding FAPAR, northness,
#     and eastness, for comparison against the full analysis_num set above.
#     Outputs suffixed _FinalSet.
# ==============================================================================
finalset_num <- analysis_num[!analysis_num %in%
                             c("fapar_annual_fapar", "aspect_northness", "aspect_eastness")]

cor_mat_fs <- cor(analysis_data[, finalset_num], method = "spearman",
                  use = "pairwise.complete.obs")

pairs_idx_fs <- which(upper.tri(cor_mat_fs), arr.ind = TRUE)
cor_pairs_fs <- data.frame(
  var1 = rownames(cor_mat_fs)[pairs_idx_fs[, 1]],
  var2 = colnames(cor_mat_fs)[pairs_idx_fs[, 2]],
  rho  = round(cor_mat_fs[pairs_idx_fs], 3)
)
cor_pairs_fs <- cor_pairs_fs[order(-abs(cor_pairs_fs$rho)), ]
flagged_pairs_fs <- cor_pairs_fs[abs(cor_pairs_fs$rho) >= COLLINEAR_THRESH, ]

cat(sprintf("\n=== [FinalSet] Collinearity (Spearman |rho| >= %.1f) ===\n", COLLINEAR_THRESH))
cat(sprintf("Flagged pairs: %d\n\n", nrow(flagged_pairs_fs)))
print(flagged_pairs_fs, row.names = FALSE)

hc_ord_fs     <- hclust(as.dist(1 - cor_mat_fs))$order
var_levels_fs <- rownames(cor_mat_fs)[hc_ord_fs]

cor_df_fs <- data.frame(
  var1 = rep(rownames(cor_mat_fs), ncol(cor_mat_fs)),
  var2 = rep(colnames(cor_mat_fs), each = nrow(cor_mat_fs)),
  rho  = as.vector(cor_mat_fs)
)
cor_df_fs$var1 <- factor(plot_labels[as.character(cor_df_fs$var1)],
                         levels = plot_labels[var_levels_fs])
cor_df_fs$var2 <- factor(plot_labels[as.character(cor_df_fs$var2)],
                         levels = plot_labels[var_levels_fs])

jpeg(file.path(outdir, "env_correlation_heatmap_FinalSet.jpg"),
     width = 5, height = 4.5, units = "in", res = 300)
print(
  ggplot(cor_df_fs, aes(x = var2, y = var1, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title   = element_blank())
)
dev.off()
cat(sprintf("[FinalSet] Heatmap written to %s\n",
            file.path(outdir, "env_correlation_heatmap_FinalSet.jpg")))

pca_data_fs <- analysis_data[, finalset_num]
pca_data_fs <- pca_data_fs[complete.cases(pca_data_fs), ]
pca_sids_fs <- analysis_data$SampleID[complete.cases(analysis_data[, finalset_num])]
cat(sprintf("\n[FinalSet] PCA: %d variables, %d complete-case samples\n",
            ncol(pca_data_fs), nrow(pca_data_fs)))

pca_res_fs <- prcomp(pca_data_fs, center = TRUE, scale. = TRUE)

var_exp_fs <- summary(pca_res_fs)$importance
cat("\n[FinalSet] Variance explained:\n")
print(round(var_exp_fs[, 1:min(10, ncol(var_exp_fs))], 3))

n_pc_fs <- min(10, ncol(pca_res_fs$rotation))
cat("\n[FinalSet] Loadings (PC1–PC", n_pc_fs, "):\n", sep = "")
print(round(pca_res_fs$rotation[, 1:n_pc_fs], 3))

pct_var_fs  <- 100 * pca_res_fs$sdev^2 / sum(pca_res_fs$sdev^2)
scree_df_fs <- data.frame(PC = seq_along(pct_var_fs), pct_var = pct_var_fs,
                          cumvar = cumsum(pct_var_fs))

jpeg(file.path(outdir, "env_pca_scree_FinalSet.jpg"), width = 900, height = 700)
print(
  ggplot(scree_df_fs, aes(x = PC, y = pct_var)) +
    geom_col(fill = "steelblue") +
    geom_line(aes(y = cumvar)) +
    geom_point(aes(y = cumvar)) +
    geom_hline(yintercept = 80, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = seq_len(nrow(scree_df_fs))) +
    labs(x = "PC", y = "Variance explained (%)") +
    theme_bw(base_size = 14)
)
dev.off()
cat(sprintf("[FinalSet] Scree plot written to %s\n",
            file.path(outdir, "env_pca_scree_FinalSet.jpg")))

scores_fs <- as.data.frame(pca_res_fs$x)
scores_fs$SampleID <- pca_sids_fs
scores_fs <- merge(scores_fs, merged[, c("SampleID", "Clusters")], by = "SampleID")

arrows_df_fs  <- as.data.frame(pca_res_fs$rotation[, 1:2])
arrows_df_fs$var <- plot_labels[rownames(arrows_df_fs)]

jpeg(file.path(outdir, "env_pca_biplot_FinalSet.jpg"), width = 1100, height = 1000)
print(
  ggplot(scores_fs, aes(x = PC1, y = PC2, colour = Clusters)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_segment(data = arrows_df_fs,
                 aes(x = 0, y = 0,
                     xend = PC1 * scale_arrow, yend = PC2 * scale_arrow),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "grey30") +
    geom_text(data = arrows_df_fs,
              aes(x = PC1 * scale_arrow * 1.1, y = PC2 * scale_arrow * 1.1,
                  label = var),
              inherit.aes = FALSE, size = 12, colour = "grey20") +
    scale_colour_manual(values = popcolors) +
    labs(x = sprintf("PC1 (%.1f%%)", pct_var_fs[1]),
         y = sprintf("PC2 (%.1f%%)", pct_var_fs[2])) +
    theme_bw(base_size = 70) +
    theme(legend.position = "none")
)
dev.off()
cat(sprintf("[FinalSet] Biplot written to %s\n",
            file.path(outdir, "env_pca_biplot_FinalSet.jpg")))
