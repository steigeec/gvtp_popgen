# historic_data.R
#
# Companion to modern_data.R, following the undergrad's histenvvarsmodels.R
# in spirit (long-window fog via Kling instead of MODIS/Werner). Uses
# three _historic environmental variables — BIO18 (summer precip), BIO15
# (precip seasonality), and Kling fog (the real observed 1950-2014 fog
# series; the analog of "fog" here, since MODIS FLCC only starts 2000) —
# each a rolling 10-year mean ending in (and including) the specimen's
# own germination year, not a fixed reference window (see
# extract_pheno_env_data.R Section 15). Kling's 1950-2014 coverage caps
# this script's usable range at germination years 1959-2014 — a small
# sample (~20 specimens) accepted per user decision 2026-09-06, rather
# than blending in a different fog source for 2015+ specimens (which
# would just be modern_data.R's fog_modern again, i.e. that script's
# job). Answers a different question than modern_data.R: does long-term
# fog/precip character predict these traits in older specimens, vs.
# modern_data.R's recent-decade climate in specimens from 2009-2024.
#
# Against three floret-count/proportion response variables (no plant
# size/shape traits — those aren't part of the undergrad's historic
# analyses):
#   Count_disk_flourets, Count_ray_flourets, ray_prop (ray / total florets)
# ranks the candidate models by AICc, and writes:
#   historic_data_output/model_selection_report.csv
#     — one block per phenotype, the best 3 models by AICc weight, with a
#       row per fixed-effect term (p-value, estimate +/- SE, z/t-score, N).
#       Carries a header note caveating that p-values are not corrected
#       for multiple comparisons — this is exploratory model selection.
#   historic_data_output/response_plot_<phenotype>_rank<N>.png
#     — for any of those top-3 models where the intercept AND every fixed
#       effect is significant (p < 0.05): a continuous raster + a few
#       labeled contour lines (color = model prediction), masked to the
#       observed data's convex hull; a single-predictor model instead gets
#       a 1D response curve; a 3+ predictor model facets across the
#       low/median/high (10th/50th/90th percentile) of its next-strongest
#       predictor (see modern_data.R for why — same logic, copied as-is).
#   historic_data_output/ray_disk_over_time.png
#     — disk count, ray count, and ray proportion vs. collection year
#       (full historic-to-modern range), each with a linear trend line
#       (raw/unadjusted — see year_trend_report.csv for the locality-
#       adjusted version).
#   historic_data_output/year_trend_report.csv
#     — trait ~ year_z + (1 | Locality) for each phenotype: does year
#       predict the trait once site is accounted for (Section 8)?
#
# Env variables are fit on z-scored (mean 0, SD 1) versions, additive plus
# pairwise interactions only (MAX_INTERACTION_ORDER below) — modern_data.R
# found higher-order interactions among these correlated env variables
# produce spurious, collinearity-driven "best" models, so the same cap
# applies here. Response-surface plots refit the winning formula on the
# raw (unscaled) variables for interpretable axes.
#
# Poisson fits (floret counts) are checked for overdispersion and
# automatically refit as negative binomial when overdispersed — see
# modern_data.R for the full rationale, identical here.
#
# Count_ray_flourets is a fusion of dissection counts (pheno_file) and
# iNat-photo counts (inatrays.csv), justified by the same paired Wilcoxon
# test as modern_data.R (re-run here since it's on the full, not
# year-filtered, phenotype file — same result either way, since the test
# itself was never year-restricted in the first place).
#
# All localities are kept, including small ones — see modern_data.R's
# locality-exclusion discussion; same reasoning applies here.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(glmmTMB)
  library(ggplot2)
  library(viridis)
  library(geometry)
  library(isoband)
  library(patchwork)  # side-by-side panel composition (Section 3)
})

# ============================================================
# Paths
# ============================================================

env_file      <- "/private/groups/meyerlab/esteiger/pheno/emma_gee_results_final.csv"
pheno_file    <- "/private/groups/meyerlab/esteiger/pheno/gtinatandvars.csv"
inatrays_file <- "/private/groups/meyerlab/esteiger/pheno/inatrays.csv"
out_dir       <- "/private/groups/meyerlab/esteiger/pheno/historic_data_output"
dir.create(out_dir, showWarnings = FALSE)

SIG_ALPHA <- 0.05
N_TOP     <- 3

AXIS_LABELS <- c(Kling_historic = "Fog", BIO18_historic = "Warmest quarter precip.",
                  BIO15_historic = "Precip. seasonality")
PLOT_BASE_SIZE     <- 20
CONTOUR_LABEL_SIZE <- 7
CONTOUR_LABEL_PT   <- CONTOUR_LABEL_SIZE * ggplot2::.pt

# ============================================================
# 1. Ray-floret count: validate + fuse dissection and iNat-photo counts
#    (identical logic/result to modern_data.R — see there for the full
#    rationale; not year-restricted in either script)
# ============================================================

pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
pheno$recordNumber_CollectorNumber <- trimws(as.character(pheno$recordNumber_CollectorNumber))
pheno$year <- suppressWarnings(as.numeric(pheno$year))
pheno$Count_ray_flourets <- suppressWarnings(as.numeric(pheno$Count_ray_flourets))
pheno$iNat_ray_count     <- suppressWarnings(as.numeric(pheno$iNat_ray_count))

paired <- pheno %>%
  transmute(Count_ray_flourets, iNat_ray_count) %>%
  drop_na(Count_ray_flourets, iNat_ray_count)
cat(sprintf("Ray-count validation: %d specimens with both dissection and iNat counts\n", nrow(paired)))
if (nrow(paired) >= 3) {
  diffs     <- paired$Count_ray_flourets - paired$iNat_ray_count
  shapiro_p <- shapiro.test(diffs)$p.value
  wt <- wilcox.test(paired$Count_ray_flourets, paired$iNat_ray_count, paired = TRUE)
  cat(sprintf("  Shapiro-Wilk on differences: p = %.2g (%s)\n", shapiro_p,
              if (shapiro_p < 0.05) "non-normal, paired Wilcoxon used" else "normal"))
  cat(sprintf("  Paired Wilcoxon signed-rank: V = %.0f, p = %.3f (mean dissection = %.2f, mean iNat = %.2f)\n",
              wt$statistic, wt$p.value, mean(paired$Count_ray_flourets), mean(paired$iNat_ray_count)))
  if (wt$p.value < SIG_ALPHA) {
    warning("Dissection and iNat ray counts differ significantly (p = ", round(wt$p.value, 3),
            ") — fusing them below is NOT well-justified by this test.")
  } else {
    cat("  Not significant — dissection and iNat ray counts are treated as interchangeable.\n\n")
  }
}

normalize_locality <- function(x) {
  x <- gsub("_", "", x)
  as.character(ifelse(is.na(x) | x == "", NA_character_, x))
}

# ============================================================
# 2. Load + merge phenotype and _historic environmental data
#
#    No explicit year filter — env-data availability does the
#    restricting instead: Kling_historic/BIO18_historic/BIO15_historic
#    are only non-NA for germination years 1959-2014 (the rolling
#    10-years-inclusive-of-germination-year window's coverage limits,
#    bound by Kling's real 1950-2014 data — see extract_pheno_env_data.R
#    Section 15), a small sample (~20 specimens) accepted per user
#    decision 2026-09-06 rather than blending in a different fog source
#    for 2015+ specimens. Join key/merge logic otherwise identical to
#    modern_data.R.
# ============================================================

env <- read.csv(env_file, stringsAsFactors = FALSE)
env$SampleID <- trimws(as.character(env$SampleID))
env <- distinct(env[, c("SampleID", "BIO18_historic", "BIO15_historic", "Kling_historic")])

historic <- pheno %>%
  transmute(
    SampleID             = recordNumber_CollectorNumber,
    year,
    Locality             = normalize_locality(Locality),
    Count_disk_flourets  = suppressWarnings(as.numeric(Count_disk_flourets)),
    Count_ray_flourets
  ) %>%
  distinct()

cat(sprintf("Phenotype records (any year): %d\n", nrow(historic)))

historic <- inner_join(historic, env, by = "SampleID")
cat(sprintf("After merging with %s: %d records (%d dropped — no matching SampleID)\n",
            basename(env_file), nrow(historic),
            sum(!pheno$recordNumber_CollectorNumber %in% env$SampleID)))

# --- iNat-only observations (no physical specimen, ray count only) ---

inat <- read.csv(inatrays_file, stringsAsFactors = FALSE)
inat$SampleID <- trimws(as.character(inat$inat_Specimen))
inat_expansion <- inat %>%
  filter(!is.na(iNat_ray_count), !is.na(latitude), !is.na(longitude),
         !(SampleID %in% historic$SampleID)) %>%
  transmute(SampleID, year = Year, Locality = normalize_locality(Locality),
            Count_disk_flourets = NA_real_, Count_ray_flourets = iNat_ray_count) %>%
  distinct() %>%
  inner_join(env, by = "SampleID")
cat(sprintf("iNat-only ray-count records added (inatrays.csv, no matching pheno_file specimen): %d\n",
            nrow(inat_expansion)))

historic <- bind_rows(historic, inat_expansion)

na_loc <- is.na(historic$Locality)
if (any(na_loc)) {
  inat_loc_by_id <- setNames(normalize_locality(inat$Locality), inat$SampleID)
  historic$Locality[na_loc] <- unname(inat_loc_by_id[historic$SampleID[na_loc]])
  cat(sprintf("Backfilled Locality for %d record(s) from inatrays.csv (still NA: %d)\n",
              sum(!is.na(historic$Locality[na_loc])), sum(is.na(historic$Locality))))
}

historic$Locality <- droplevels(as.factor(historic$Locality))
historic$ray_prop <- with(historic, Count_ray_flourets / (Count_disk_flourets + Count_ray_flourets))
historic$Count_total_florets <- with(historic, Count_disk_flourets + Count_ray_flourets)

# z-scored predictors used for model selection
env_vars <- c(Fog = "Kling_historic", BIO18 = "BIO18_historic", BIO15 = "BIO15_historic")
for (nm in names(env_vars)) {
  historic[[paste0(nm, "_z")]] <- as.numeric(scale(historic[[env_vars[[nm]]]]))
}
z_vars <- paste0(names(env_vars), "_z")   # c("Fog_z","BIO18_z","BIO15_z")

# ============================================================
# 3. Ray:disk floret change over time
#
#    Uses the full phenotype file regardless of _historic env-data
#    availability (unlike everything else below, which is restricted to
#    germination years 1959-2014) — this is a simple temporal-trend
#    check of collection year itself, requested alongside the env-driven
#    analysis below.
# ============================================================

time_data <- historic %>%
  drop_na(year) %>%
  transmute(year, Disk = Count_disk_flourets, Ray = Count_ray_flourets) %>%
  pivot_longer(-year, names_to = "metric", values_to = "value") %>%
  drop_na(value) %>%
  mutate(metric = factor(metric, levels = c("Disk", "Ray")))

# Sized up from the standard CONTOUR_LABEL_PT — this plot is meant to read
# at a smaller display scale (e.g. inset/thumbnail), so its embedded font
# needs to be oversized relative to the figure's own (narrower) dimensions.
RD_TIME_FONT_PT <- CONTOUR_LABEL_PT * 1.3

# --- BIO18-over-time panel (left half of the combined figure below),
#     replicating plot_env_time.R's "coord_overall_loess" version — a
#     single trimmed LOESS fit pooling every individual coordinate's
#     annual BIO18 value (env_wi_range.csv, full 1927-2024 PRISM record),
#     not broken out by locality — themed (no title/gridlines, enlarged
#     fonts, tightened legend spacing) and colored (black line, grey
#     ribbon) to match the disk/ray panel's own single-line look. Trim
#     left at the original script's span/2-tied fraction (BIO18's
#     span=0.5, chosen there for extra smoothing against BIO18's
#     zero-inflated spikes) rather than extended toward the present —
#     pushing the trim out to reach 2024 made the line visibly wigglier
#     without a corresponding gain in what it's meant to show here.
env_range <- read.csv("/private/groups/meyerlab/esteiger/pheno/env_wi_range.csv", stringsAsFactors = FALSE)

bio18_long <- env_range %>%
  select(Locality, Lat, Long, matches("^BIO18_[0-9]{4}$")) %>%
  pivot_longer(matches("^BIO18_[0-9]{4}$"), names_to = "year", values_to = "value") %>%
  mutate(year = as.integer(sub("BIO18_", "", year))) %>%
  drop_na(value)

# Trimmed LOESS: at each end, a fit draws on data from only one side of the
# window, so the outer span/2 fraction of the year range is dropped before
# plotting (identical logic to plot_env_time.R's trim_loess()).
trim_loess_bio18 <- function(year, value, span = 0.5) {
  df <- data.frame(year = year, value = value)
  if (nrow(df) < 5) return(data.frame())
  yr_range <- range(df$year)
  trim <- (span / 2) * diff(yr_range)
  rel_min <- yr_range[1] + trim; rel_max <- yr_range[2] - trim
  if (rel_min >= rel_max) return(data.frame())
  fit <- tryCatch(loess(value ~ year, data = df, span = span), error = function(e) NULL)
  if (is.null(fit)) return(data.frame())
  pred_yr <- seq(yr_range[1], yr_range[2], by = 1)
  pred <- suppressWarnings(predict(fit, newdata = data.frame(year = pred_yr), se = TRUE))
  data.frame(year = pred_yr, fitted = pred$fit,
             ymin = pred$fit - 1.96 * pred$se.fit, ymax = pred$fit + 1.96 * pred$se.fit) %>%
    filter(year >= rel_min & year <= rel_max & !is.na(fitted))
}

bio18_loess <- trim_loess_bio18(bio18_long$year, bio18_long$value)

p_bio18 <- ggplot(bio18_loess, aes(x = year, y = fitted)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4) +
  geom_line(color = "black", linewidth = 0.6) +
  labs(x = "Year", y = "Warmest quarter precip.") +
  theme_bw(base_size = PLOT_BASE_SIZE) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = RD_TIME_FONT_PT),
        axis.text  = element_text(size = RD_TIME_FONT_PT * 0.75),
        axis.title.x = element_text(margin = margin(t = 4)),
        plot.margin = margin(t = 5, r = 8, b = 2, l = 5))

p_time <- ggplot(time_data, aes(x = year, y = value, shape = metric, color = metric)) +
  # small jitter only, so many specimens sharing an exact (year, count) pair
  # (common in recent, densely-sampled years) don't fully overlap into one
  # dot — geom_smooth() below is unaffected, since jitter is applied only
  # to this layer's drawn positions, not to the underlying data it fits on.
  geom_point(position = position_jitter(width = 0.3, height = 0.3, seed = 42),
             alpha = 0.85, size = 2.3) +
  # color is a fixed override here (not mapped) so both trend lines render
  # in the same neutral black/grey — only the points are colored per metric
  geom_smooth(method = "lm", color = "black", fill = "grey70", linewidth = 0.6, se = TRUE) +
  scale_shape_manual(values = c("Disk" = 16, "Ray" = 17)) +
  scale_color_manual(values = c("Disk" = "grey65", "Ray" = "black")) +
  labs(x = "Year", y = "Floret count", shape = NULL, color = NULL) +
  theme_bw(base_size = PLOT_BASE_SIZE) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = RD_TIME_FONT_PT),
        axis.text  = element_text(size = RD_TIME_FONT_PT * 0.75),
        axis.title.x = element_text(margin = margin(t = 4)),
        legend.text = element_text(size = RD_TIME_FONT_PT * 0.75),
        legend.position = "bottom",
        legend.box.spacing = unit(2, "pt"),
        legend.margin = margin(t = 0, b = 0),
        plot.margin = margin(t = 5, r = 8, b = 2, l = 5))
ggsave(file.path(out_dir, "ray_disk_over_time.png"), p_time, width = 6.5, height = 7, dpi = 150)
cat(sprintf("\nRay/disk-over-time plot written to %s\n", file.path(out_dir, "ray_disk_over_time.png")))

# Captured for the combined a/b/c/d figure at the very end of this script
# (panel_c = BIO18 LOESS, panel_d = ray/disk-over-time) — carries forward
# PHENOTYPE_DATA_PANELS$a/$b already set in the modern section above.
if (!exists("PHENOTYPE_DATA_PANELS")) PHENOTYPE_DATA_PANELS <- list()
PHENOTYPE_DATA_PANELS$c <- p_bio18
PHENOTYPE_DATA_PANELS$d <- p_time

p_combined <- p_bio18 + p_time
ggsave(file.path(out_dir, "env_and_florets_over_time.png"), p_combined, width = 13, height = 7, dpi = 150)
cat(sprintf("Combined BIO18 + ray/disk-over-time figure written to %s\n",
            file.path(out_dir, "env_and_florets_over_time.png")))

# ============================================================
# 4. Generic helpers: build formula set, fit, AICc-rank, extract coefs
#    (identical to modern_data.R — see there for the full rationale)
# ============================================================

MAX_INTERACTION_ORDER <- 2

build_formulas <- function(response, vars) {
  # Null model included explicitly (2026-09-06, per user request) so the
  # AICc ranking directly shows whether ANY env-variable combination
  # actually outperforms "no predictors at all" — relevant here since
  # historic_data.R's sample is small enough that this isn't a given.
  null_model <- paste(response, "~ 1 + (1 | Locality)")
  additive <- unlist(lapply(seq_along(vars), function(i) {
    combn(vars, i, function(v) paste(response, "~", paste(v, collapse = " + "), "+ (1 | Locality)"))
  }))
  max_order <- min(MAX_INTERACTION_ORDER, length(vars))
  interactive <- unlist(lapply(2:max_order, function(i) {
    combn(vars, i, function(v) paste(response, "~", paste(v, collapse = " * "), "+ (1 | Locality)"))
  }))
  unique(c(null_model, additive, interactive))
}

dispersion_stat <- function(model) {
  rp <- residuals(model, type = "pearson")
  sum(rp^2) / df.residual(model)
}
DISPERSION_THRESH <- 1.4

fit_one <- function(formula_str, data, family) {
  suppressWarnings(suppressMessages(tryCatch({
    m <- switch(family,
      poisson  = glmer(as.formula(formula_str), data = data, family = poisson,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))),
      gamma    = glmer(as.formula(formula_str), data = data, family = Gamma(link = "log"),
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))),
      gaussian = lmerTest::lmer(as.formula(formula_str), data = data),
      beta     = glmmTMB(as.formula(formula_str), data = data, family = beta_family(link = "logit"))
    )
    if (identical(family, "poisson")) {
      disp <- dispersion_stat(m)
      if (is.finite(disp) && disp > DISPERSION_THRESH) {
        cat(sprintf("    overdispersed (%.2f > %.1f), refitting as negative binomial: %s\n",
                    disp, DISPERSION_THRESH, formula_str))
        m <- glmmTMB(as.formula(formula_str), data = data, family = nbinom2(link = "log"))
      }
    }
    m
  }, error = function(e) NULL)))
}

fit_all <- function(response, data, family, vars = z_vars) {
  formulas <- build_formulas(response, vars)
  models <- lapply(formulas, fit_one, data = data, family = family)
  ok <- !sapply(models, is.null)
  vcov_finite <- function(m) {
    v <- vcov(m)
    if (is.list(v)) v <- v$cond
    all(is.finite(as.matrix(v)))
  }
  ok[ok] <- sapply(models[ok], vcov_finite)
  cat(sprintf("  %s: %d / %d candidate models fit successfully\n",
              response, sum(ok), length(formulas)))
  list(models = models[ok], formulas = formulas[ok])
}

aicc_rank <- function(models) {
  aic <- sapply(models, AIC)
  k   <- sapply(models, function(m) attr(logLik(m), "df"))
  n   <- sapply(models, function(m) nobs(m))
  aicc  <- aic + (2 * k * (k + 1)) / (n - k - 1)
  delta <- aicc - min(aicc)
  w     <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))
  data.frame(index = seq_along(models), K = k, AICc = aicc, DeltaAICc = delta, Weight = w) %>%
    arrange(AICc)
}

get_coefs <- function(model) {
  if (inherits(model, "glmmTMB")) {
    cc <- as.data.frame(summary(model)$coefficients$cond)
    names(cc) <- c("estimate", "std_error", "statistic", "p_value")
  } else if (inherits(model, "lmerModLmerTest")) {
    cc <- as.data.frame(summary(model)$coefficients)
    names(cc) <- c("estimate", "std_error", "df", "statistic", "p_value")
    cc$df <- NULL
  } else {
    cc <- as.data.frame(summary(model)$coefficients)
    names(cc) <- c("estimate", "std_error", "statistic", "p_value")
  }
  cc$term <- rownames(cc)
  rownames(cc) <- NULL
  cc[, c("term", "estimate", "std_error", "statistic", "p_value")]
}

display_name <- function(term) {
  term %>%
    str_replace_all(setNames(names(env_vars), z_vars)) %>%
    str_replace_all(":", " x ")
}

# ============================================================
# 5. Response-surface / response-curve plotting
#    (identical to modern_data.R — see there for the full rationale)
# ============================================================

mask_to_hull <- function(grid_xy, obs_xy) {
  obs_xy <- obs_xy[complete.cases(obs_xy), ]
  hull <- convhulln(as.matrix(obs_xy))
  inhulln(hull, as.matrix(grid_xy))
}

compute_contour_lines <- function(x_seq, y_seq, z, n_levels = 5) {
  if (all(is.na(z))) return(NULL)
  breaks <- pretty(range(z, na.rm = TRUE), n = n_levels)
  breaks <- breaks[breaks > min(z, na.rm = TRUE) & breaks < max(z, na.rm = TRUE)]
  if (length(breaks) == 0) return(NULL)
  zmat <- t(matrix(z, nrow = length(x_seq), ncol = length(y_seq)))
  iso  <- isoband::isolines(x_seq, y_seq, zmat, levels = breaks)
  lines_df <- bind_rows(lapply(names(iso), function(lvl) {
    seg <- iso[[lvl]]
    if (length(seg$x) == 0) return(NULL)
    data.frame(x = seg$x, y = seg$y, id = paste(lvl, seg$id), level = as.numeric(lvl))
  }))
  if (is.null(lines_df) || nrow(lines_df) == 0) return(NULL)
  labels_df <- lines_df %>%
    group_by(level) %>%
    slice(which.min(abs(x - median(x)))) %>%
    ungroup()
  list(lines = lines_df, labels = labels_df)
}

plot_response_surface <- function(response, formula_str, data, family, phenotype_label, out_path) {
  raw_formula_str <- formula_str
  for (i in seq_along(z_vars)) raw_formula_str <- gsub(z_vars[i], env_vars[[i]], raw_formula_str, fixed = TRUE)

  raw_model <- fit_one(raw_formula_str, data, family)
  if (is.null(raw_model)) {
    warning("Raw-scale refit failed for: ", raw_formula_str, " — skipping plot.")
    return(invisible(NULL))
  }

  predictors <- all.vars(as.formula(raw_formula_str))
  predictors <- predictors[predictors %in% unname(env_vars)]

  plot_label <- if (family == "gaussian") sub("\\(log ", "(", phenotype_label) else phenotype_label

  pred_fun <- function(newdata) {
    if (inherits(raw_model, "glmmTMB")) {
      predict(raw_model, newdata = newdata, type = "response", re.form = NA)
    } else if (family == "gaussian") {
      exp(predict(raw_model, newdata = newdata, re.form = NA))
    } else {
      predict(raw_model, newdata = newdata, re.form = NA, type = "response")
    }
  }

  if (length(predictors) == 1) {
    xv <- predictors[1]
    rng <- range(data[[xv]], na.rm = TRUE)
    grid <- data.frame(x = seq(rng[1], rng[2], length.out = 200))
    names(grid) <- xv
    grid$pred <- pred_fun(grid)

    p <- ggplot(grid, aes(x = .data[[xv]], y = pred)) +
      geom_line(linewidth = 1.2, color = "#238A8DFF") +
      geom_rug(data = data, aes(x = .data[[xv]]), sides = "b", alpha = 0.4, inherit.aes = FALSE) +
      labs(x = AXIS_LABELS[[xv]], y = paste("Predicted", plot_label)) +
      theme_bw(base_size = PLOT_BASE_SIZE) +
      theme(panel.grid = element_blank(),
            axis.title = element_text(size = CONTOUR_LABEL_PT),
            axis.text  = element_text(size = CONTOUR_LABEL_PT))
    ggsave(out_path, p, width = 7, height = 5.5, dpi = 150)
    return(invisible(p))
  }

  coefs <- get_coefs(raw_model) %>% filter(term %in% predictors)
  axis_vars <- coefs$term[order(-abs(coefs$statistic))][1:2]
  held_vars <- setdiff(predictors, axis_vars)

  rng1 <- range(data[[axis_vars[1]]], na.rm = TRUE)
  rng2 <- range(data[[axis_vars[2]]], na.rm = TRUE)
  base_grid <- expand.grid(
    setNames(list(seq(rng1[1], rng1[2], length.out = 100),
                  seq(rng2[1], rng2[2], length.out = 100)),
             axis_vars)
  )
  base_grid$in_hull <- as.vector(mask_to_hull(base_grid[, axis_vars], data[, axis_vars]))

  if (length(held_vars) == 0) {
    grid <- base_grid
    grid$pred <- pred_fun(grid)
    grid$pred_masked <- ifelse(grid$in_hull, grid$pred, NA)
  } else {
    held_coefs   <- coefs %>% filter(term %in% held_vars)
    facet_var    <- held_coefs$term[order(-abs(held_coefs$statistic))][1]
    mean_held    <- setdiff(held_vars, facet_var)
    for (hv in mean_held) base_grid[[hv]] <- mean(data[[hv]], na.rm = TRUE)

    facet_vals   <- quantile(data[[facet_var]], c(0.1, 0.5, 0.9), na.rm = TRUE)
    facet_var_label <- AXIS_LABELS[[facet_var]]
    facet_labels <- c(sprintf("%s: %.3g (10th pct)", facet_var_label, facet_vals[1]),
                       sprintf("%s: %.3g (median)",   facet_var_label, facet_vals[2]),
                       sprintf("%s: %.3g (90th pct)", facet_var_label, facet_vals[3]))

    grid <- bind_rows(lapply(seq_along(facet_vals), function(i) {
      g <- base_grid
      g[[facet_var]]   <- facet_vals[[i]]
      g$facet_label    <- factor(facet_labels[i], levels = facet_labels)
      g
    }))
    grid$pred <- pred_fun(grid)
    grid$pred_masked <- ifelse(grid$in_hull, grid$pred, NA)
  }

  x_seq <- sort(unique(base_grid[[axis_vars[1]]]))
  y_seq <- sort(unique(base_grid[[axis_vars[2]]]))
  if (length(held_vars) == 0) {
    co <- compute_contour_lines(x_seq, y_seq, grid$pred_masked)
    contour_lines  <- co$lines
    contour_labels <- co$labels
  } else {
    facet_groups <- levels(grid$facet_label)
    co_list <- lapply(facet_groups, function(fl) {
      sub <- grid[grid$facet_label == fl, ]
      co  <- compute_contour_lines(x_seq, y_seq, sub$pred_masked)
      if (is.null(co)) return(NULL)
      co$lines$facet_label  <- factor(fl, levels = facet_groups)
      co$labels$facet_label <- factor(fl, levels = facet_groups)
      co
    })
    contour_lines  <- bind_rows(lapply(co_list, `[[`, "lines"))
    contour_labels <- bind_rows(lapply(co_list, `[[`, "labels"))
  }

  p <- ggplot() +
    geom_raster(data = grid, aes(x = .data[[axis_vars[1]]], y = .data[[axis_vars[2]]], fill = pred_masked))
  if (!is.null(contour_lines) && nrow(contour_lines) > 0) {
    p <- p + geom_path(data = contour_lines, aes(x = x, y = y, group = id),
                        color = "white", linewidth = 0.4, alpha = 0.8)
  }
  if (!is.null(contour_labels) && nrow(contour_labels) > 0) {
    p <- p + geom_label(data = contour_labels, aes(x = x, y = y, label = signif(level, 3)),
                         size = CONTOUR_LABEL_SIZE, linewidth = 0, fill = "white", alpha = 0.75,
                         color = "grey20", fontface = "bold")
  }
  p <- p +
    geom_point(data = data, aes(x = .data[[axis_vars[1]]], y = .data[[axis_vars[2]]]), size = 2.5) +
    scale_fill_viridis_c(guide = "none", na.value = "transparent") +
    coord_cartesian(expand = FALSE) +
    labs(x = AXIS_LABELS[[axis_vars[1]]], y = AXIS_LABELS[[axis_vars[2]]]) +
    theme_bw(base_size = PLOT_BASE_SIZE) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = CONTOUR_LABEL_PT),
          axis.text  = element_text(size = CONTOUR_LABEL_PT))
  if (length(held_vars) > 0) p <- p + facet_wrap(~facet_label)
  out_width <- if (length(held_vars) > 0) 15 else 8
  ggsave(out_path, p, width = out_width, height = 6.5, dpi = 150)
  invisible(p)
}

# ============================================================
# 6-7. Per-phenotype env-variable GLM (cross-sectional: does trait LEVEL
#      relate to env LEVEL across specimens) + its report — RE-ENABLED
#      2026-09-06 per user decision, while Section 10's mediation test is
#      temporarily disabled (see below) during work on modern_data.R's
#      parallel Werner-fog mediation section. Re-disable by wrapping in
#      if (FALSE) { ... } again if the mediation focus resumes.
# ============================================================

if (TRUE) {

report_rows <- list()

run_phenotype <- function(response, phenotype_label, family, data) {
  cat(sprintf("\n=== %s (%s) ===\n", phenotype_label, response))
  # NOT applying modern_data.R's uniform-N-across-candidates fix here.
  # z_vars includes Fog (Kling_historic), whose real coverage is ~12-16
  # specimens (vs. ~159-162 for BIO18/BIO15) — forcing every candidate
  # model in the set to share one N would crush BIO18/BIO15-only models
  # down to Kling's tiny sample too, for no benefit (tested 2026-09-07:
  # collapsed ray_prop's whole candidate set to n=12). Per user decision,
  # historic keeps its original mixed-N behavior — each formula's N is
  # whatever glmer/glmmTMB's own per-formula NA handling gives it — and
  # accepts the resulting cross-model AICc caveat as the lesser cost.
  data$Locality <- droplevels(data$Locality)
  fit <- fit_all(response, data, family)
  if (length(fit$models) == 0) {
    warning("No candidate models converged for ", response, " — skipping.")
    return(invisible(NULL))
  }
  ranked <- aicc_rank(fit$models) %>% slice_head(n = N_TOP)

  for (rank in seq_len(nrow(ranked))) {
    idx   <- ranked$index[rank]
    model <- fit$models[[idx]]
    fstr  <- fit$formulas[[idx]]
    coefs <- get_coefs(model)

    intercept_ok <- coefs$p_value[coefs$term == "(Intercept)"] < SIG_ALPHA
    var_coefs    <- coefs %>% filter(term != "(Intercept)")
    all_sig      <- intercept_ok && length(var_coefs$p_value) > 0 &&
                     all(!is.na(var_coefs$p_value)) && all(var_coefs$p_value < SIG_ALPHA)

    formula_display <- if (grepl("~ 1 \\+ \\(1", fstr)) {
      "(null model)"
    } else {
      paste(sapply(strsplit(gsub(".*~ (.*) \\+ \\(1.*", "\\1", fstr), " ")[[1]],
                   function(w) if (w %in% c("+", "*")) w else display_name(w)),
            collapse = " ")
    }

    all_terms <- bind_rows(coefs[coefs$term == "(Intercept)", ], var_coefs)
    n_obs <- nobs(model)
    for (r in seq_len(nrow(all_terms))) {
      report_rows[[length(report_rows) + 1]] <<- data.frame(
        Phenotype    = phenotype_label,
        N            = n_obs,
        Formula      = if (r == 1) formula_display else "",
        Model_weight = if (r == 1) round(ranked$Weight[rank], 3) else NA,
        AICc         = if (r == 1) round(ranked$AICc[rank], 2) else NA,
        Delta_AICc   = if (r == 1) round(ranked$DeltaAICc[rank], 2) else NA,
        Variable     = if (all_terms$term[r] == "(Intercept)") "(Intercept)" else display_name(all_terms$term[r]),
        p_value      = signif(all_terms$p_value[r], 2),
        Estimate_SE  = sprintf("%.3f ± %.3f", all_terms$estimate[r], all_terms$std_error[r]),
        statistic    = round(all_terms$statistic[r], 3)
      )
    }

    if (all_sig) {
      cat(sprintf("  Rank %d model is fully significant — plotting.\n", rank))
      out_path <- file.path(out_dir, sprintf("response_plot_%s_rank%d.png", response, rank))
      plot_response_surface(response, fstr, data, family, phenotype_label, out_path)
    }
  }
}

run_phenotype("Count_disk_flourets", "Count disk florets", "poisson",
              historic %>% drop_na(Count_disk_flourets))
run_phenotype("Count_ray_flourets", "Count ray florets", "poisson",
              historic %>% drop_na(Count_ray_flourets))
run_phenotype("Count_total_florets", "Count total florets", "poisson",
              historic %>% drop_na(Count_total_florets))
run_phenotype("ray_prop", "Proportion ray florets", "beta",
              historic %>% drop_na(ray_prop))

# ============================================================
# 7. Write report
# ============================================================

report <- bind_rows(report_rows) %>%
  rename(`Model weight` = Model_weight, `AICc` = AICc, `Delta AICc` = Delta_AICc,
         `p-value` = p_value, `Estimate (β) ± std error` = Estimate_SE, `z/t-score` = statistic)

report_path <- file.path(out_dir, "model_selection_report.csv")

n_candidate_models <- length(build_formulas("y", z_vars))
caveat <- sprintf(
  "# Exploratory model selection: %d candidate models (all additive + pairwise-interaction combinations of %s) fit per phenotype, ranked by AICc. p-values are NOT corrected for multiple comparisons across models or terms — treat as hypothesis-generating, not confirmatory.",
  n_candidate_models, paste(names(env_vars), collapse = "/"))
writeLines(caveat, report_path)
suppressWarnings(write.table(report, report_path, sep = ",", row.names = FALSE, na = "",
                              append = TRUE, qmethod = "double"))
cat(sprintf("\nReport written to %s (%d rows, plus a multiple-comparisons caveat header)\n",
            report_path, nrow(report)))

} # end if (TRUE) — Sections 6-7 re-enabled

# ============================================================
# 8. Formal test: does each trait change with collection year, after
#    accounting for locality?
#
#    The ray_disk_over_time.png plot (Section 3) is a raw/unadjusted
#    trend — it doesn't distinguish "traits genuinely change over time"
#    from "different eras happened to sample different sites". This
#    fits trait ~ year_z + (1 | Locality) for each phenotype (same
#    family/dispersion-check/fusion logic as the env-variable models
#    above) to test the year effect with site pooled out.
# ============================================================

cat("\n=== Year-trend models (trait ~ year, locality as random intercept) ===\n\n")

year_trend_rows <- list()
run_year_trend <- function(response, phenotype_label, family, data) {
  data <- data %>% drop_na(year)
  data$Locality <- droplevels(data$Locality)
  data$year_z <- as.numeric(scale(data$year))
  m <- fit_one(paste(response, "~ year_z + (1 | Locality)"), data, family)
  if (is.null(m)) {
    warning("Year-trend model failed to converge for ", response)
    return(invisible(NULL))
  }
  cc <- get_coefs(m) %>% filter(term == "year_z")
  cat(sprintf("  %s (n=%d): year β = %.4f ± %.4f, z/t = %.3f, p = %.3g%s\n",
              phenotype_label, nobs(m), cc$estimate, cc$std_error, cc$statistic, cc$p_value,
              if (cc$p_value < SIG_ALPHA) "  *" else ""))
  year_trend_rows[[length(year_trend_rows) + 1]] <<- data.frame(
    Phenotype = phenotype_label, N = nobs(m), `Estimate (β) ± std error` = sprintf("%.4f ± %.4f", cc$estimate, cc$std_error),
    `z/t-score` = round(cc$statistic, 3), `p-value` = signif(cc$p_value, 2), check.names = FALSE)
}
run_year_trend("Count_disk_flourets", "Count disk florets", "poisson", historic %>% drop_na(Count_disk_flourets))
run_year_trend("Count_ray_flourets", "Count ray florets", "poisson", historic %>% drop_na(Count_ray_flourets))
run_year_trend("Count_total_florets", "Count total florets", "poisson", historic %>% drop_na(Count_total_florets))
run_year_trend("ray_prop", "Proportion ray florets", "beta", historic %>% drop_na(ray_prop))

year_trend_report <- bind_rows(year_trend_rows)
year_trend_path <- file.path(out_dir, "year_trend_report.csv")
write.csv(year_trend_report, year_trend_path, row.names = FALSE)
cat(sprintf("\nYear-trend report written to %s\n", year_trend_path))

# ============================================================
# 9. Do the env variables themselves trend over time?
#
#    A different question from the GLM above (which asks whether trait
#    LEVEL relates to env LEVEL across specimens) and from Section 8
#    (whether the trait trends with collection year) — this asks whether
#    fog/BIO18/BIO15 themselves have a secular trend, using every year of
#    data at every sampled coordinate (env_wi_range.csv: per-coordinate
#    annual values, BIO12/15/18 1927-2024, fog 2000-2025, kling_fog
#    1950-2014) rather than each specimen's own rolling 10yr-inclusive
#    window — far more power for THIS question, since it isn't
#    bottlenecked by which specimens happen to have a computable window.
#    Locality-capped (15/coordinate, matching the project's standard
#    collinearity-check methodology) since Strauss has 27 unique
#    coordinates here vs. Dangermond's 2.
#
#    If an env variable shows no year trend, it can't be mediating the
#    year effects found in Section 8. If it does, that's necessary but
#    not sufficient for a mediating role — genuine mediation would need
#    trait ~ year + env to show a smaller/less-significant year term
#    than trait ~ year alone (Section 8); not attempted here given how
#    few specimens have both a trait value and a fully-covered env
#    window to test that on.
# ============================================================

cat("\n=== Env-variable-over-time trend (env_wi_range.csv, locality-capped) ===\n\n")

range_env <- read.csv("/private/groups/meyerlab/esteiger/pheno/env_wi_range.csv",
                       stringsAsFactors = FALSE, check.names = FALSE)

set.seed(42)
range_capped <- do.call(rbind, lapply(split(range_env, range_env$Locality), function(d) {
  if (nrow(d) > 15) d[sample(nrow(d), 15), ] else d
}))
cat(sprintf("Locality-capped coordinates: %d (from %d total)\n", nrow(range_capped), nrow(range_env)))

pivot_env_long <- function(data, var_prefix) {
  cols <- grep(paste0("^", var_prefix, "_[0-9]{4}$"), names(data), value = TRUE)
  yrs  <- as.integer(sub(paste0("^", var_prefix, "_"), "", cols))
  out <- data.frame(Locality = rep(data$Locality, length(cols)),
                     Lat = rep(data$Lat, length(cols)), Long = rep(data$Long, length(cols)),
                     year = rep(yrs, each = nrow(data)),
                     value = unlist(data[, cols], use.names = FALSE))
  out
}

env_trend_rows <- list()
run_env_trend <- function(var_prefix, var_label) {
  long <- pivot_env_long(range_capped, var_prefix) %>% drop_na(value)
  long$Locality <- droplevels(as.factor(long$Locality))
  long$year_z <- as.numeric(scale(long$year))
  m <- lmerTest::lmer(value ~ year_z + (1 | Locality), data = long)
  cc <- get_coefs(m) %>% filter(term == "year_z")
  cat(sprintf("  %s (n=%d obs, %d coords x %d yrs %d-%d): year beta = %.4f +/- %.4f, t = %.3f, p = %.3g%s\n",
              var_label, nobs(m), length(unique(long$Lat)), length(unique(long$year)),
              min(long$year), max(long$year), cc$estimate, cc$std_error, cc$statistic, cc$p_value,
              if (cc$p_value < SIG_ALPHA) "  *" else ""))
  env_trend_rows[[length(env_trend_rows) + 1]] <<- data.frame(
    Variable = var_label, N_obs = nobs(m), N_coords = length(unique(paste(long$Lat, long$Long))),
    `Estimate (β) ± std error` = sprintf("%.4f ± %.4f", cc$estimate, cc$std_error),
    `t-score` = round(cc$statistic, 3), `p-value` = signif(cc$p_value, 2), check.names = FALSE)
}
run_env_trend("kling_fog", "Fog (Kling)")
run_env_trend("BIO18", "BIO18 (summer precip)")
run_env_trend("BIO15", "BIO15 (precip seasonality)")

env_trend_report <- bind_rows(env_trend_rows)
env_trend_path <- file.path(out_dir, "env_trend_report.csv")
write.csv(env_trend_report, env_trend_path, row.names = FALSE)
cat(sprintf("\nEnv-trend report written to %s\n", env_trend_path))

# ============================================================
# 10. Mediation test: do Fog and BIO18's trends over time explain the
#     year-trend in ray floret proportion? — TEMPORARILY DISABLED
#     2026-09-06 while building the parallel Werner-fog mediation section
#     in modern_data.R; Sections 6-7 (the cross-sectional GLM) are
#     re-enabled above in the meantime. Re-enable by wrapping this
#     section back in if (TRUE) { ... } once that work is done.
#
#     modern_data.R found Fog and BIO18 (summer precip) both positively
#     predict ray_prop (n=138, Fog beta=+0.052 p=0.036, BIO18 beta=+0.064
#     p=0.041). Section 9 above found Fog is DECLINING over time
#     (p=9e-25) while BIO18 is INCREASING (p=1.8e-5) — so naively, a
#     declining Fog would predict a declining ray_prop, yet Section 8
#     found ray_prop trending upward (marginally, p=0.056). BIO18's rise
#     is the candidate explanation: if its (positive) effect on ray_prop
#     outweighs Fog's declining contribution, the net trend could still
#     be upward. Tested here by comparing the year coefficient in
#     ray_prop ~ year alone vs. ray_prop ~ year + Fog + BIO18 (does it
#     attenuate once the env variables are controlled for?), plus a
#     BIO18-only version (more power — BIO18_historic's PRISM-based
#     coverage is broader than Fog_historic's Kling-based one, so more
#     specimens qualify when Fog isn't required).
#
#     All using _historic (20yr rolling, inclusive of germination year)
#     variables, restricted to whichever specimens have them — a small,
#     likely underpowered sample; treat as suggestive, not confirmatory.
# ============================================================

if (FALSE) {

run_mediation <- function(response, label, family) {
  cat(sprintf("\n--- %s ---\n\n", label))

  med_full_data <- historic %>% drop_na(!!sym(response), year, Fog_z, BIO18_z)
  med_full_data$Locality <- droplevels(med_full_data$Locality)
  med_full_data$year_z <- as.numeric(scale(med_full_data$year))
  cat(sprintf("Full model (Fog + BIO18): n = %d specimens\n", nrow(med_full_data)))

  m_year_only_full <- fit_one(paste(response, "~ year_z + (1 | Locality)"), med_full_data, family)
  m_full           <- fit_one(paste(response, "~ year_z + Fog_z + BIO18_z + (1 | Locality)"), med_full_data, family)

  if (!is.null(m_year_only_full) && !is.null(m_full)) {
    cc_yo <- get_coefs(m_year_only_full) %>% filter(term == "year_z")
    cat(sprintf("  year alone (same n=%d, for fair comparison): beta=%.4f se=%.4f z=%.3f p=%.3g\n",
                nobs(m_year_only_full), cc_yo$estimate, cc_yo$std_error, cc_yo$statistic, cc_yo$p_value))
    cc_full <- get_coefs(m_full)
    for (r in seq_len(nrow(cc_full))) {
      cat(sprintf("  full model %-12s beta=%.4f se=%.4f z=%.3f p=%.3g%s\n",
                  paste0(cc_full$term[r], ":"), cc_full$estimate[r], cc_full$std_error[r],
                  cc_full$statistic[r], cc_full$p_value[r],
                  if (cc_full$p_value[r] < SIG_ALPHA) "  *" else ""))
    }
  } else {
    cat("  Full model (year + Fog + BIO18) failed to converge on this small sample.\n")
  }

  # BIO18-only version — more specimens (doesn't need Fog/Kling), directly
  # tests the "BIO18's rise alone" hypothesis
  med_bio18_data <- historic %>% drop_na(!!sym(response), year, BIO18_z)
  med_bio18_data$Locality <- droplevels(med_bio18_data$Locality)
  med_bio18_data$year_z <- as.numeric(scale(med_bio18_data$year))
  cat(sprintf("\nBIO18-only model: n = %d specimens\n", nrow(med_bio18_data)))

  m_year_only_b18 <- fit_one(paste(response, "~ year_z + (1 | Locality)"), med_bio18_data, family)
  m_bio18         <- fit_one(paste(response, "~ year_z + BIO18_z + (1 | Locality)"), med_bio18_data, family)

  if (!is.null(m_year_only_b18) && !is.null(m_bio18)) {
    cc_yo_b <- get_coefs(m_year_only_b18) %>% filter(term == "year_z")
    cat(sprintf("  year alone (same n=%d, for fair comparison): beta=%.4f se=%.4f z=%.3f p=%.3g\n",
                nobs(m_year_only_b18), cc_yo_b$estimate, cc_yo_b$std_error, cc_yo_b$statistic, cc_yo_b$p_value))
    cc_b18 <- get_coefs(m_bio18)
    for (r in seq_len(nrow(cc_b18))) {
      cat(sprintf("  +BIO18 model %-12s beta=%.4f se=%.4f z=%.3f p=%.3g%s\n",
                  paste0(cc_b18$term[r], ":"), cc_b18$estimate[r], cc_b18$std_error[r],
                  cc_b18$statistic[r], cc_b18$p_value[r],
                  if (cc_b18$p_value[r] < SIG_ALPHA) "  *" else ""))
    }
  } else {
    cat("  BIO18-only model failed to converge.\n")
  }
  invisible(NULL)
}

cat("\n=== Section 10: Mediation test — do Fog/BIO18's trends over time explain each phenotype's year-trend? ===\n")

run_mediation("ray_prop", "Proportion ray florets", "beta")
run_mediation("Count_ray_flourets", "Count ray florets", "poisson")
run_mediation("Count_disk_flourets", "Count disk florets", "poisson")
run_mediation("Count_total_florets", "Count total florets", "poisson")

} # end if (FALSE) — Section 10 (mediation) temporarily disabled

# ============================================================
# 11. Combined a/b/c/d figure — single row:
#       a: Height:width ratio response curve (modern, rank3 top significant model)
#       b: Proportion ray florets response surface (modern, rank1 top significant model)
#       c: BIO18-over-time, pooled LOESS (historic, env_wi_range.csv)
#       d: Ray/disk floret counts over time (historic)
#     Panels a/b default to CONTOUR_LABEL_PT sizing (tuned for their own
#     standalone response_plot_*.png files); overridden here to match c/d's
#     larger RD_TIME_FONT_PT so all four panels share one font size in
#     this combined figure specifically — their individually-saved files
#     are untouched, since this re-themes captured copies, not the
#     originals.
# ============================================================

combined_out_dir <- "/private/groups/meyerlab/esteiger/pheno/phenotype_data_output"
dir.create(combined_out_dir, showWarnings = FALSE)

if (all(c("a", "b", "c", "d") %in% names(PHENOTYPE_DATA_PANELS))) {
  match_font <- function(p) {
    p + theme(axis.title = element_text(size = RD_TIME_FONT_PT),
              axis.text  = element_text(size = RD_TIME_FONT_PT * 0.75))
  }
  panel_a <- match_font(PHENOTYPE_DATA_PANELS$a)
  panel_b <- match_font(PHENOTYPE_DATA_PANELS$b)
  panel_c <- match_font(PHENOTYPE_DATA_PANELS$c)
  panel_d <- match_font(PHENOTYPE_DATA_PANELS$d)

  combined_4panel <- (panel_a | panel_b | panel_c | panel_d) +
    plot_annotation(tag_levels = "a")
  ggsave(file.path(combined_out_dir, "combined_4panel_a_b_c_d.png"), combined_4panel,
         width = 26, height = 6.5, dpi = 150)
  cat(sprintf("\nCombined 4-panel figure written to %s\n",
              file.path(combined_out_dir, "combined_4panel_a_b_c_d.png")))
} else {
  warning("Not all 4 panels (a/b/c/d) were captured — combined figure skipped. Missing: ",
          paste(setdiff(c("a", "b", "c", "d"), names(PHENOTYPE_DATA_PANELS)), collapse = ", "))
}
