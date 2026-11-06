#!/usr/bin/env Rscript
.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")

# ---- Replicate data pipeline from inbreeding_plot.R -------------------------
inbreeding <- read.table("ngsF_results.txt", header=T)
inbreeding$SampleID <- gsub("_rmdup.bam", "", inbreeding$BAM)
pop <- read.table("/private/groups/meyerlab/esteiger/field_data.txt", header=T, sep="\t")
pop$SampleID <- as.character(pop$SampleID)
inbreeding <- merge(inbreeding, pop, by="SampleID")
inbreeding$log <- log10(inbreeding$Cov)
inbreeding <- inbreeding[!(inbreeding$SampleID %in% c("4531")), ]
inbmed <- inbreeding[inbreeding$log > 0.85, ]
inbmed$Year <- as.factor(inbmed$Year)

mod_cov <- lm(NGSF ~ log10(Cov), data=inbmed)
inbmed$residual <- residuals(mod_cov)

# ---- Merge env data (Year 2020, good coords only) ---------------------------
env <- read.table("/private/groups/meyerlab/esteiger/env_var.txt", header=T, sep="\t")
inbmed$SampleID <- trimws(as.character(inbmed$SampleID))
env$SampleID    <- trimws(as.character(env$SampleID))
env <- merge(inbmed, env, by="SampleID")
env <- env[env$Year.x == "2020", ]

edna_env <- read.csv("/private/groups/meyerlab/esteiger/eDNA_env_data.csv", header=TRUE)
edna_env$SampleID <- trimws(as.character(edna_env$SampleID))
env <- merge(env, edna_env, by="SampleID", suffixes=c("", "_edna"))

# ---- Site-level aggregation -------------------------------------------------
site <- aggregate(residual ~ Clusters, data=env, FUN=mean)

Ne_vals <- c(CDFW=42, Dangermond=165, Point_Conception=827,
             Strauss=712, VAFB_Arguello=165, VAFB_Lion_Head=455)
site$Ne <- Ne_vals[as.character(site$Clusters)]

# ---- Drought variable list --------------------------------------------------
drought_vars <- c(
  # from env_var.txt
  "Percent_volume_soil_water_content_at_0_cm",
  "PercVolSoilH20_10cm",
  "BIO12", "BIO15", "BIO18",
  # from eDNA (soil water)
  "soil_soil_water_content_v_percent_0_30cm",
  "soil_soil_water_content_v_percent_0cm",
  # precipitation / chirps
  "chirps_heavy_precip_days_count_30y",
  "chirps_rainfall_intensity_sdii_mm_day_30y",
  "chirps_heavy_precip_days_count_10y",
  "chirps_rainfall_intensity_sdii_mm_day_10y",
  # terraclimate 30y
  "terraclimate_climate_water_deficit_mm_30y",
  "terraclimate_palmer_drought_severity_index_30y",
  "terraclimate_precipitation_accumulation_mm_30y",
  "terraclimate_snow_water_equivalent_mm_30y",
  "terraclimate_soil_moisture_mm_30y",
  "terraclimate_vapor_pressure_deficit_kpa_30y",
  "terraclimate_vapor_pressure_kpa_30y",
  "terraclimate_annual_precip_30y",
  "terraclimate_precip_driest_month_30y",
  "terraclimate_precip_driest_quarter_30y",
  "terraclimate_precip_wettest_month_30y",
  "terraclimate_precip_wettest_quarter_30y",
  "terraclimate_precip_seasonality_30y",
  "terraclimate_precip_annual_cv_30y",
  "terraclimate_driest_year_precip_mm_30y",
  "terraclimate_severe_drought_month_count_30y",
  "terraclimate_aridity_index_30y",
  "terraclimate_soil_moisture_seasonality_cv_30y",
  # terraclimate 10y
  "terraclimate_climate_water_deficit_mm_10y",
  "terraclimate_palmer_drought_severity_index_10y",
  "terraclimate_precipitation_accumulation_mm_10y",
  "terraclimate_snow_water_equivalent_mm_10y",
  "terraclimate_soil_moisture_mm_10y",
  "terraclimate_vapor_pressure_deficit_kpa_10y",
  "terraclimate_vapor_pressure_kpa_10y",
  "terraclimate_annual_precip_10y",
  "terraclimate_precip_driest_month_10y",
  "terraclimate_precip_driest_quarter_10y",
  "terraclimate_precip_wettest_month_10y",
  "terraclimate_precip_wettest_quarter_10y",
  "terraclimate_precip_seasonality_10y",
  "terraclimate_precip_annual_cv_10y",
  "terraclimate_driest_year_precip_mm_10y",
  "terraclimate_severe_drought_month_count_10y",
  "terraclimate_aridity_index_10y",
  "terraclimate_soil_moisture_seasonality_cv_10y"
)

# Keep only variables that exist in merged env
drought_vars <- drought_vars[drought_vars %in% names(env)]
cat("Testing", length(drought_vars), "variables\n\n")

# ---- Aggregate each drought var to site level, then test --------------------
results <- data.frame(
  variable       = character(),
  r_inb          = numeric(),   # simple Pearson r vs inb residual
  p_inb          = numeric(),
  r_Ne           = numeric(),   # Pearson r vs Ne
  p_Ne           = numeric(),
  r_partial      = numeric(),   # added-variable r (controlling for Ne)
  p_partial      = numeric(),
  stringsAsFactors = FALSE
)

for (v in drought_vars) {
  vals <- env[[v]]
  # skip if non-numeric or all-NA or all-same
  if (!is.numeric(vals)) next
  if (sum(!is.na(vals)) < 3) next

  # aggregate to site mean
  site_v <- aggregate(env[[v]] ~ env$Clusters, FUN=mean, na.rm=TRUE)
  names(site_v) <- c("Clusters", "dvar")
  d <- merge(site, site_v, by="Clusters")
  if (nrow(d) < 4 || sd(d$dvar, na.rm=TRUE) == 0) next

  # simple correlation with inb residual
  ct_inb  <- tryCatch(cor.test(d$dvar, d$residual, method="pearson"), error=function(e) NULL)
  # correlation with Ne
  ct_Ne   <- tryCatch(cor.test(d$dvar, d$Ne,       method="pearson"), error=function(e) NULL)
  # added-variable: partial r after controlling for Ne
  e_inb_Ne  <- tryCatch(residuals(lm(residual ~ Ne, data=d)), error=function(e) NULL)
  e_dvar_Ne <- tryCatch(residuals(lm(dvar     ~ Ne, data=d)), error=function(e) NULL)
  ct_part   <- tryCatch(cor.test(e_dvar_Ne, e_inb_Ne, method="pearson"), error=function(e) NULL)

  if (is.null(ct_inb) || is.null(ct_Ne) || is.null(ct_part)) next

  results <- rbind(results, data.frame(
    variable  = v,
    r_inb     = round(ct_inb$estimate,  3),
    p_inb     = round(ct_inb$p.value,   4),
    r_Ne      = round(ct_Ne$estimate,   3),
    p_Ne      = round(ct_Ne$p.value,    4),
    r_partial = round(ct_part$estimate, 3),
    p_partial = round(ct_part$p.value,  4),
    stringsAsFactors = FALSE
  ))
}

# ---- Sort by |r_partial| and print all, highlight top 3 --------------------
results <- results[order(-abs(results$r_partial)), ]
rownames(results) <- NULL

cat("=================================================================\n")
cat("ALL RESULTS (sorted by |partial r| controlling for Ne)\n")
cat("=================================================================\n")
print(results, row.names=FALSE)

cat("\n=================================================================\n")
cat("TOP 3 PERFORMERS (by |partial r| after removing Ne effect)\n")
cat("=================================================================\n")

for (i in 1:3) {
  v   <- results$variable[i]
  d_v <- aggregate(env[[v]] ~ env$Clusters, FUN=mean, na.rm=TRUE)
  names(d_v) <- c("Clusters", "dvar")
  d   <- merge(site, d_v, by="Clusters")

  e_inb_Ne  <- residuals(lm(residual ~ Ne, data=d))
  e_dvar_Ne <- residuals(lm(dvar     ~ Ne, data=d))
  m_full    <- lm(residual ~ dvar + Ne, data=d)
  m_simple  <- lm(residual ~ dvar,      data=d)

  cat(sprintf("\n--- #%d: %s ---\n", i, v))
  cat(sprintf("  Simple:  r = %.3f, p = %.4f\n",    results$r_inb[i], results$p_inb[i]))
  cat(sprintf("  vs Ne:   r = %.3f, p = %.4f\n",    results$r_Ne[i],  results$p_Ne[i]))
  cat(sprintf("  Partial: r = %.3f, p = %.4f\n",    results$r_partial[i], results$p_partial[i]))
  cat("  lm(residual ~ dvar):\n")
  print(summary(m_simple)$coefficients)
  cat(sprintf("  R² = %.3f\n", summary(m_simple)$r.squared))
  cat("  lm(residual ~ dvar + Ne):\n")
  print(summary(m_full)$coefficients)
  cat(sprintf("  R² = %.3f\n\n", summary(m_full)$r.squared))
}
