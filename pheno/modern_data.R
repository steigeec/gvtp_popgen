# modern_data.R
#
# Reworked, self-contained version of the undergrad's modern-specimen GLM
# analyses (see gvtp_popgen/pheno/modern_models.R, heightmodels.R,
# updatedvarmodels.R for the original, ad hoc versions this replaces).
#
# For every phenotype record with a computable env window (effectively
# germination year 2009-2024 — see Section 2), fits all additive and
# interaction combinations of two _modern environmental variables — fog
# and BIO18 (summer precip), each a rolling 10-year mean ending in (and
# including) the specimen's own germination year, not a fixed reference
# window (see extract_pheno_env_data.R Section 15). BIO15 (precip
# seasonality) was dropped 2026-09-06 — collinear with fog (Spearman rho
# -0.84 to -0.89 depending on window length) beyond this project's 0.7
# threshold — against seven response variables:
#   Count_disk_flourets, Count_ray_flourets, Count_total_florets (disk + ray),
#   ray_prop (ray / total florets — the pistillate/fertile fraction; disk
#   florets are functionally staminate), Height_in, Width_in, Height_width
#   (height:width ratio)
# ranks the candidate models by AICc, and writes:
#   modern_data_output/model_selection_report.csv
#     — one block per phenotype, the best 3 models by AICc weight, with a
#       row per fixed-effect term (p-value, estimate +/- SE, z/t-score).
#       Carries a header note caveating that p-values are not corrected
#       for multiple comparisons (many candidate models/terms) — this is
#       exploratory model selection, not a confirmatory analysis.
#   modern_data_output/response_plot_<phenotype>_<rank>.png
#     — for any of those top-3 models where the intercept AND every fixed
#       effect is significant (p < 0.05): a 2D contour plot (color = model
#       prediction) over its two strongest predictors, masked to the
#       observed data's convex hull; a single-predictor model instead gets
#       a 1D response curve; a 3+ predictor model facets across the
#       low/median/high (10th/50th/90th percentile) of its next-strongest
#       predictor rather than collapsing it to one mean-held slice (any
#       remaining held predictors beyond that are held at their mean).
#
# Env variables are fit on z-scored (mean 0, SD 1) versions for model
# comparability and to keep interaction terms numerically well-behaved;
# response-surface plots refit the winning formula on the raw (unscaled)
# variables purely so the plot axes are in natural, interpretable units
# (mirrors the original modern_models.R workflow: select on scaled
# predictors, plot on raw ones).
#
# Poisson fits (floret counts) are checked for overdispersion (Pearson
# dispersion statistic > 1.4, the same threshold the undergrad's scripts
# used manually via blmeco::dispersion_glmer) and automatically refit as
# negative binomial (glmmTMB nbinom2) when overdispersed.
#
# Count_ray_flourets is a fusion of dissection counts (pheno_file) and
# iNat-photo counts (inatrays.csv, for specimens with no physical
# collection — ray count only, never disk count), justified by a paired
# Wilcoxon test on double-measured specimens (see Section 1, run first)
# showing no significant difference between the two methods.
#
# All localities are kept, including small ones (Dangermond n=4, CDFW n=6,
# Agujas n=7, CoastalLompoc n=1 — a single specimen backfilled from
# inatrays.csv, see Section 2). The undergrad's original scripts dropped
# Dangermond/CDFW in most models for statistical reasons (small n); that's
# exactly what (1 | Locality) as a random effect is designed to handle —
# a small group just contributes less information and shrinks toward the
# grand mean, rather than needing to be excluded (confirmed empirically:
# every phenotype's full 21-model candidate set still converges with all
# 8 localities included, including the n=1 singleton).

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lme4)
  library(lmerTest)   # Satterthwaite p-values for the Gaussian (lmer) models
  library(glmmTMB)
  library(ggplot2)
  library(viridis)
  library(geometry)   # convex hull masking for contour plots
  library(isoband)    # contour-line computation for the labeled raster plots
})

# ============================================================
# Paths
# ============================================================

env_file      <- "/private/groups/meyerlab/esteiger/pheno/emma_gee_results_final.csv"
pheno_file    <- "/private/groups/meyerlab/esteiger/pheno/gtinatandvars.csv"
inatrays_file <- "/private/groups/meyerlab/esteiger/pheno/inatrays.csv"
out_dir       <- "/private/groups/meyerlab/esteiger/pheno/modern_data_output"
dir.create(out_dir, showWarnings = FALSE)

SIG_ALPHA <- 0.05
N_TOP     <- 3

# Plot styling: no legend (contour-line labels already give exact values),
# no title/subtitle, no gridlines, human-readable axis names, larger text
# throughout (contour-line labels bigger still, since they carry the
# entire "what number is this" job with the legend gone).
AXIS_LABELS <- c(fog_modern = "Fog", BIO18_modern = "Warmest quarter precip.")
PLOT_BASE_SIZE     <- 20
CONTOUR_LABEL_SIZE <- 7
# geom_label()'s `size` is in mm; convert to the pt size it actually
# renders at (via ggplot2's own mm-to-pt factor) so axis title/text can be
# set to match it exactly, rather than relying on theme_bw()'s default
# (smaller, relative-to-base-size) axis text scaling.
CONTOUR_LABEL_PT <- CONTOUR_LABEL_SIZE * ggplot2::.pt

# ============================================================
# 1. Ray-floret count: validate + fuse dissection and iNat-photo counts
#
#     The undergrad's stattests.R set up (but didn't capture output for)
#     a paired comparison of dissection-based Count_ray_flourets vs.
#     iNat_ray_count for specimens measured both ways, to justify using
#     one as a stand-in for the other. Re-run here — on the full
#     phenotype file, not just modern (year > 2000) records, matching
#     the scope their own script used — as a paired Wilcoxon signed-rank
#     test (Shapiro-Wilk on the differences rejects normality, matching
#     their own "#not normal" comment, so a paired t-test isn't
#     appropriate): if not significant, dissection and iNat ray counts
#     are treated as interchangeable for the rest of this script — use
#     dissection where available, fall back to the iNat count otherwise.
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

# inatrays.csv spells locality names with underscores (e.g. "Point_Conception",
# "VAFB_Lion_head"); pheno_file doesn't (e.g. "PointConception",
# "VAFBLionHead") — normalize both so the same physical site collapses to
# one random-effect level instead of silently splitting into two. Also
# turns "" into NA (one inatrays.csv row has no locality recorded).
normalize_locality <- function(x) {
  x <- gsub("_", "", x)
  # as.character() guards against base ifelse()'s zero-length quirk: when
  # `test` has length 0, ifelse() returns logical(0) regardless of the
  # yes/no argument types, which breaks bind_rows() when inat_expansion
  # legitimately has 0 rows (now common, since Section 13 closes most of
  # the env-data gap upstream, leaving little for this fallback to add).
  as.character(ifelse(is.na(x) | x == "", NA_character_, x))
}

# ============================================================
# 2. Load + merge phenotype and _modern environmental data
#
#    pheno_file (gtinatandvars.csv) is the undergrad's original phenotype
#    source file (received directly from them 2026-09-06; previously this
#    script used updatedvarsandpheno.csv, a merged bridging table someone
#    had built from this same file + an older emma_gee_results_final.csv
#    pass — same 224 specimens, same columns, just without 2 duplicate
#    rows and without the now-superseded env columns). Its own BIO12/15/18/
#    fog/CHILI columns predate the current _germyr/_25yr/_<N>yr convention,
#    so they're ignored here — only phenotype measurements are used; env
#    values are joined fresh from the current emma_gee_results_final.csv.
#
#    Join key: emma_gee_results_final.csv$SampleID ==
#              gtinatandvars.csv$recordNumber_CollectorNumber
#    (both are the herbarium/collector record number, or the bare iNat
#    observation ID for iNat-only records; confirmed by matching against
#    all_coord_gvtp.csv's Dataset == "PHENO" rows). Not every phenotype
#    record matches — a few have a collector number recorded as a date
#    (excluded, not a real ID) or reference a coordinate the extraction
#    pipeline hadn't processed yet. That coordinate gap (~50 records,
#    including several the undergrad's original scripts manually excluded
#    as "weird" — e.g. "233000618", the 8503-8505 range — which turned
#    out to just be missing env data, not arbitrary drops) is closed by
#    extract_pheno_env_data.R Section 13; re-run that script first if
#    emma_gee_results_final.csv predates it.
# ============================================================

env <- read.csv(env_file, stringsAsFactors = FALSE)
env$SampleID <- trimws(as.character(env$SampleID))
env <- distinct(env[, c("SampleID", "fog_modern", "BIO18_modern")])

# No explicit year filter here anymore — env-data availability (fog_modern/
# BIO18_modern are only non-NA for germYear 2009-2024, the rolling
# 10-years-inclusive-of-germination-year window's coverage limits;
# see extract_pheno_env_data.R Section 15) naturally restricts the sample
# to that range, superseding the old fixed "year > 2000" rule.
modern <- pheno %>%
  filter(!is.na(year)) %>%
  transmute(
    SampleID             = recordNumber_CollectorNumber,
    year,
    Locality             = normalize_locality(Locality),
    Count_disk_flourets  = suppressWarnings(as.numeric(Count_disk_flourets)),
    Count_ray_flourets,
    Height_in            = suppressWarnings(as.numeric(Height_in)),
    Width_in             = suppressWarnings(as.numeric(Width_in)),
    Height_width         = suppressWarnings(as.numeric(Height_width))
  ) %>%
  distinct()

cat(sprintf("Modern (year > 2000) phenotype records: %d\n", nrow(modern)))

is_modern <- !is.na(pheno$year) & pheno$year > 2000
modern <- inner_join(modern, env, by = "SampleID")
cat(sprintf("After merging with %s: %d records (%d dropped — no matching SampleID)\n",
            basename(env_file), nrow(modern),
            sum(!pheno$recordNumber_CollectorNumber[is_modern] %in% env$SampleID)))

# --- iNat-only observations (inatrays.csv — no physical specimen exists,
#     so only a ray count is possible, never a disk-floret count), folded
#     in as additional Count_ray_flourets-only records, excluding any
#     SampleID already represented above to avoid double-counting a
#     specimen that was measured both ways ---

inat <- read.csv(inatrays_file, stringsAsFactors = FALSE)
inat$SampleID <- trimws(as.character(inat$inat_Specimen))
inat_expansion <- inat %>%
  filter(!is.na(iNat_ray_count), !is.na(Year), Year > 2000, !is.na(latitude), !is.na(longitude),
         !(SampleID %in% modern$SampleID)) %>%
  transmute(SampleID, year = Year, Locality = normalize_locality(Locality), Count_disk_flourets = NA_real_,
            Count_ray_flourets = iNat_ray_count,
            Height_in = NA_real_, Width_in = NA_real_, Height_width = NA_real_) %>%
  distinct() %>%
  inner_join(env, by = "SampleID")
cat(sprintf("iNat-only ray-count records added (inatrays.csv, no matching pheno_file specimen): %d\n",
            nrow(inat_expansion)))

modern <- bind_rows(modern, inat_expansion)

# Backfill a still-missing Locality (e.g. specimen 233000618, blank in
# pheno_file) from inatrays.csv when available, rather than losing that
# row to an NA random-effect group — small localities are kept, not
# dropped, elsewhere in this script (see the locality-exclusion
# discussion), so a singleton shouldn't be discarded either.
na_loc <- is.na(modern$Locality)
if (any(na_loc)) {
  inat_loc_by_id <- setNames(normalize_locality(inat$Locality), inat$SampleID)
  modern$Locality[na_loc] <- unname(inat_loc_by_id[modern$SampleID[na_loc]])
  cat(sprintf("Backfilled Locality for %d record(s) from inatrays.csv (still NA: %d)\n",
              sum(!is.na(modern$Locality[na_loc])), sum(is.na(modern$Locality))))
}

modern$Locality <- droplevels(as.factor(modern$Locality))
# Ray florets are pistillate/fertile (seed-producing), disk florets are
# functionally staminate — ray_prop is the biologically meaningful
# "fertile floret" fraction, and the framing this project's pollinator-
# attraction/sex-allocation literature comparisons use. Mathematically
# ray_prop = 1 - disk_prop, so this is a sign-flip of the same beta-
# regression fit, not a different result (logit(ray_prop) = -logit(disk_prop)).
modern$ray_prop <- with(modern, Count_ray_flourets / (Count_disk_flourets + Count_ray_flourets))
modern$Count_total_florets <- with(modern, Count_disk_flourets + Count_ray_flourets)

# z-scored predictors used for model selection. BIO15 was dropped (was
# BIO15_modern) 2026-09-06: Fog x BIO15 collinearity (Spearman rho -0.84
# at a 10yr rolling window, -0.89 at 20yr — tested longer windows
# specifically to see if it would resolve this, it didn't) exceeded the
# project's 0.7 locality-capped threshold regardless of window length.
env_vars <- c(Fog = "fog_modern", BIO18 = "BIO18_modern")
for (nm in names(env_vars)) {
  modern[[paste0(nm, "_z")]] <- as.numeric(scale(modern[[env_vars[[nm]]]]))
}
z_vars <- paste0(names(env_vars), "_z")   # c("Fog_z","BIO18_z")

# ============================================================
# 3. Generic helpers: build formula set, fit, AICc-rank, extract coefs
# ============================================================

# All additive (+) combinations of `vars`, plus pairwise (*) interactions,
# each with a (1 | Locality) random intercept. The original modern_models.R
# / heightmodels.R / updatedvarmodels.R combn() loops went all the way to a
# full-order interaction; with the (formerly 4, now 3) env variables here
# notably correlated and n ~115-200, higher-order interaction terms were
# found to be practically unidentifiable (with the earlier 4-variable
# _25yr/CWD set: the full Fog*BIO18*BIO15*CWD model came out "best" by
# AICc weight ~1.0 for both height and width, with coefficient SEs 5-10x
# the estimates — collinearity-driven noise, not signal). Capped at 2-way
# interactions here for that reason; with only 3 variables now, this cap
# just excludes the single 3-way interaction term.
MAX_INTERACTION_ORDER <- 2

build_formulas <- function(response, vars) {
  additive <- unlist(lapply(seq_along(vars), function(i) {
    combn(vars, i, function(v) paste(response, "~", paste(v, collapse = " + "), "+ (1 | Locality)"))
  }))
  max_order <- min(MAX_INTERACTION_ORDER, length(vars))
  interactive <- unlist(lapply(2:max_order, function(i) {
    combn(vars, i, function(v) paste(response, "~", paste(v, collapse = " * "), "+ (1 | Locality)"))
  }))
  unique(c(additive, interactive))
}

# Pearson dispersion statistic (sum of squared Pearson residuals / residual
# df) — the same diagnostic the undergrad's modern_models.R computed
# manually per winning model (via blmeco::dispersion_glmer, "not overdispersed"
# threshold ~1.4); here it's checked automatically for every Poisson fit.
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
        # Overdispersed for a Poisson mean-variance assumption — refit as
        # negative binomial instead of trusting Poisson SEs/p-values.
        cat(sprintf("    overdispersed (%.2f > %.1f), refitting as negative binomial: %s\n",
                    disp, DISPERSION_THRESH, formula_str))
        m <- glmmTMB(as.formula(formula_str), data = data, family = nbinom2(link = "log"))
      }
    }
    m
  }, error = function(e) NULL)))
}

# Fit every formula in the candidate set; drop fits that errored or failed
# to converge to a proper covariance estimate (vcov has non-finite entries).
fit_all <- function(response, data, family, vars = z_vars) {
  formulas <- build_formulas(response, vars)
  models <- lapply(formulas, fit_one, data = data, family = family)
  ok <- !sapply(models, is.null)
  vcov_finite <- function(m) {
    v <- vcov(m)
    if (is.list(v)) v <- v$cond   # glmmTMB returns a list of component vcovs
    all(is.finite(as.matrix(v)))
  }
  ok[ok] <- sapply(models[ok], vcov_finite)
  cat(sprintf("  %s: %d / %d candidate models fit successfully\n",
              response, sum(ok), length(formulas)))
  list(models = models[ok], formulas = formulas[ok])
}

# Custom AICc table (rather than AICcmodavg::aictab, whose S3 methods don't
# reliably cover glmmTMB) built directly from the AIC()/logLik()/nobs()
# generics that lme4 and glmmTMB both implement, so counts/glm and beta
# regressions rank on the same footing.
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

# Fixed-effect coefficient table, normalized across glmer/lmerTest/glmmTMB.
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

# Display names: "fog_z" -> "Fog", "BIO18_z:CWD_z" -> "BIO18 x CWD"
display_name <- function(term) {
  term %>%
    str_replace_all(setNames(names(env_vars), z_vars)) %>%
    str_replace_all(":", " x ")
}

# ============================================================
# 4. Response-surface / response-curve plotting
#
#    Only for models where the intercept AND every fixed effect are
#    significant at SIG_ALPHA. Refits the model formula on the RAW
#    (unscaled) env variables (same data subset) so plot axes are in
#    natural units; the AICc ranking/report above is unaffected by this
#    refit, which exists only to make the plot interpretable.
# ============================================================

mask_to_hull <- function(grid_xy, obs_xy) {
  obs_xy <- obs_xy[complete.cases(obs_xy), ]
  hull <- convhulln(as.matrix(obs_xy))
  inhulln(hull, as.matrix(grid_xy))
}

# A handful of labeled contour lines over a continuous raster, computed
# directly via isoband (the same contouring engine geom_contour() uses
# internally) rather than metR/directlabels (neither installed) — avoids
# geom_contour_filled()'s stepped legend, which doesn't compress to a
# small plot (one row per bin, easily 8-10 rows). z must be in the
# expand.grid(x=x_seq, y=y_seq) row order (x fastest); returns
# list(lines=, labels=) for geom_path()/geom_label(), or NULL if z has no
# variation to contour (e.g. an entirely-masked facet).
compute_contour_lines <- function(x_seq, y_seq, z, n_levels = 5) {
  if (all(is.na(z))) return(NULL)
  breaks <- pretty(range(z, na.rm = TRUE), n = n_levels)
  breaks <- breaks[breaks > min(z, na.rm = TRUE) & breaks < max(z, na.rm = TRUE)]
  if (length(breaks) == 0) return(NULL)
  # isoband wants z as dim(length(y_seq), length(x_seq)) — i.e. z[i,j] ~
  # (x_seq[j], y_seq[i]) — the transpose of the column-major reshape that
  # matches expand.grid(x=x_seq, y=y_seq)'s row order.
  zmat <- t(matrix(z, nrow = length(x_seq), ncol = length(y_seq)))
  iso  <- isoband::isolines(x_seq, y_seq, zmat, levels = breaks)
  lines_df <- bind_rows(lapply(names(iso), function(lvl) {
    seg <- iso[[lvl]]
    if (length(seg$x) == 0) return(NULL)
    data.frame(x = seg$x, y = seg$y, id = paste(lvl, seg$id), level = as.numeric(lvl))
  }))
  if (is.null(lines_df) || nrow(lines_df) == 0) return(NULL)
  # one label per level, placed at the point closest to that level's
  # median x — a simple, reasonably central anchor without an extra
  # dependency for proper label placement
  labels_df <- lines_df %>%
    group_by(level) %>%
    slice(which.min(abs(x - median(x)))) %>%
    ungroup()
  list(lines = lines_df, labels = labels_df)
}

plot_response_surface <- function(response, formula_str, data, family, phenotype_label, out_path) {
  # rebuild the formula on raw variable names, same structure (+/interaction), same random effect
  raw_formula_str <- formula_str
  for (i in seq_along(z_vars)) raw_formula_str <- gsub(z_vars[i], env_vars[[i]], raw_formula_str, fixed = TRUE)

  raw_model <- fit_one(raw_formula_str, data, family)
  if (is.null(raw_model)) {
    warning("Raw-scale refit failed for: ", raw_formula_str, " — skipping plot.")
    return(invisible(NULL))
  }

  predictors <- all.vars(as.formula(raw_formula_str))
  predictors <- predictors[predictors %in% unname(env_vars)]

  # Gaussian models are fit on log(response) but predictions are
  # back-transformed to the raw response scale below — reflect that in the
  # plot label so it doesn't still say "(log ...)" next to raw-scale values.
  plot_label <- if (family == "gaussian") sub("\\(log ", "(", phenotype_label) else phenotype_label

  pred_fun <- function(newdata) {
    if (inherits(raw_model, "glmmTMB")) {
      predict(raw_model, newdata = newdata, type = "response", re.form = NA)
    } else if (family == "gaussian") {
      exp(predict(raw_model, newdata = newdata, re.form = NA))  # back-transform log(response)
    } else {
      predict(raw_model, newdata = newdata, re.form = NA, type = "response")
    }
  }

  if (length(predictors) == 1) {
    # --- 1D response curve ---
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

  # --- 2D contour: strongest two predictors by |z/t statistic| ---
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
    # Exactly 2 predictors — single panel, no held variables.
    grid <- base_grid
    grid$pred <- pred_fun(grid)
    grid$pred_masked <- ifelse(grid$in_hull, grid$pred, NA)
  } else {
    # 3+ predictors: facet across low/median/high (10th/50th/90th
    # percentile) of the single strongest held variable, rather than
    # collapsing it to one mean-held slice — a slice can hide an
    # interaction between an axis variable and the held one. Any further
    # held variables beyond that (rare — MAX_INTERACTION_ORDER keeps most
    # candidate models to <=3 total predictors) fall back to mean-holding
    # so the panel count stays at 3 regardless of how many are held.
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

  # Continuous raster + a few labeled contour lines, computed per facet
  # (each facet is its own prediction surface) — see compute_contour_lines().
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

# Alternative to plot_response_surface() for 2+ predictor models, added
# alongside it (not replacing it) so the two styles can be compared: an
# "interaction line plot" — the strongest predictor on the x-axis, 2-3
# colored lines for the second predictor's 10th/50th/90th percentile
# (a compact, fixed-size 3-swatch legend regardless of plot size, unlike
# a color gradient), faceted by a 3rd+ held variable exactly like
# plot_response_surface(). For single-predictor models this would just
# reproduce the same curve plot_response_surface() already made, so it's
# skipped there.
plot_response_lines <- function(response, formula_str, data, family, phenotype_label, out_path) {
  raw_formula_str <- formula_str
  for (i in seq_along(z_vars)) raw_formula_str <- gsub(z_vars[i], env_vars[[i]], raw_formula_str, fixed = TRUE)

  raw_model <- fit_one(raw_formula_str, data, family)
  if (is.null(raw_model)) {
    warning("Raw-scale refit failed for: ", raw_formula_str, " — skipping line plot.")
    return(invisible(NULL))
  }

  predictors <- all.vars(as.formula(raw_formula_str))
  predictors <- predictors[predictors %in% unname(env_vars)]
  if (length(predictors) == 1) return(invisible(NULL))  # identical to the 1D curve above

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

  coefs  <- get_coefs(raw_model) %>% filter(term %in% predictors)
  xvar   <- coefs$term[order(-abs(coefs$statistic))][1]
  colvar <- coefs$term[order(-abs(coefs$statistic))][2]
  held_vars <- setdiff(predictors, c(xvar, colvar))

  x_seq    <- seq(min(data[[xvar]], na.rm = TRUE), max(data[[xvar]], na.rm = TRUE), length.out = 200)
  col_vals <- quantile(data[[colvar]], c(0.1, 0.5, 0.9), na.rm = TRUE)
  col_labels <- sprintf("%s: %.3g (%s)", colvar, col_vals, c("10th pct", "median", "90th pct"))

  base_grid <- expand.grid(x = x_seq, col_i = seq_along(col_vals))
  names(base_grid)[1] <- xvar
  base_grid[[colvar]] <- col_vals[base_grid$col_i]
  base_grid$col_label <- factor(col_labels[base_grid$col_i], levels = col_labels)
  base_grid$col_i <- NULL
  base_grid$in_hull <- as.vector(mask_to_hull(base_grid[, c(xvar, colvar)], data[, c(xvar, colvar)]))

  if (length(held_vars) == 0) {
    grid <- base_grid
    subtitle <- NULL
  } else {
    held_coefs <- coefs %>% filter(term %in% held_vars)
    facet_var  <- held_coefs$term[order(-abs(held_coefs$statistic))][1]
    mean_held  <- setdiff(held_vars, facet_var)
    for (hv in mean_held) base_grid[[hv]] <- mean(data[[hv]], na.rm = TRUE)

    facet_vals   <- quantile(data[[facet_var]], c(0.1, 0.5, 0.9), na.rm = TRUE)
    facet_labels <- c(sprintf("%s: %.3g (10th pct)", facet_var, facet_vals[1]),
                       sprintf("%s: %.3g (median)",   facet_var, facet_vals[2]),
                       sprintf("%s: %.3g (90th pct)", facet_var, facet_vals[3]))
    grid <- bind_rows(lapply(seq_along(facet_vals), function(i) {
      g <- base_grid
      g[[facet_var]] <- facet_vals[[i]]
      g$facet_label  <- factor(facet_labels[i], levels = facet_labels)
      g
    }))
    subtitle <- if (length(mean_held) > 0) {
      sprintf("%s held at sample mean", paste(mean_held, collapse = ", "))
    } else NULL
  }

  grid$pred        <- pred_fun(grid)
  grid$pred_masked <- ifelse(grid$in_hull, grid$pred, NA)

  p <- ggplot(grid, aes(x = .data[[xvar]], y = pred_masked, color = col_label)) +
    geom_line(linewidth = 1, na.rm = TRUE) +
    scale_color_viridis_d(name = NULL, end = 0.85) +
    labs(title = sprintf("%s ~ %s", plot_label, paste(predictors, collapse = " * ")),
         subtitle = subtitle, x = xvar, y = paste("Predicted", plot_label)) +
    theme_bw(base_size = 14) + theme(legend.position = "top")
  if (length(held_vars) > 0) p <- p + facet_wrap(~facet_label)
  out_width <- if (length(held_vars) > 0) 15 else 7.5
  ggsave(out_path, p, width = out_width, height = 5.5, dpi = 150)
  invisible(p)
}

# ============================================================
# 5. Per-phenotype driver: fit candidate set, rank, report rows, plots
# ============================================================

report_rows <- list()

run_phenotype <- function(response, phenotype_label, family, data) {
  cat(sprintf("\n=== %s (%s) ===\n", phenotype_label, response))
  # Require every candidate variable (not just the ones a given formula
  # uses) to be non-missing before fitting anything, so every model in
  # the candidate set is compared on the identical set of rows. Without
  # this, a variable-poor formula (e.g. "BIO18" alone) silently retains
  # more rows than a variable-rich one (e.g. "Fog + BIO18") and can win
  # the AICc comparison on sample size rather than fit — found 2026-09-07
  # in ray_prop, where an n=159 BIO18-only model (BIO18 itself p=0.084,
  # not significant) out-weighted the n=138 Fog+BIO18 model where both
  # terms were significant, simply because it had more data to work with.
  data <- data %>% drop_na(all_of(z_vars))
  data$Locality <- droplevels(data$Locality)
  fit <- fit_all(response, data, family)
  if (length(fit$models) == 0) {
    warning("No candidate models converged for ", response, " — skipping.")
    return(invisible(NULL))
  }
  ranked <- aicc_rank(fit$models) %>% slice_head(n = N_TOP)
  plotted <- list()  # rank -> ggplot object, for any rank that gets plotted below

  for (rank in seq_len(nrow(ranked))) {
    idx   <- ranked$index[rank]
    model <- fit$models[[idx]]
    fstr  <- fit$formulas[[idx]]
    coefs <- get_coefs(model)

    intercept_ok <- coefs$p_value[coefs$term == "(Intercept)"] < SIG_ALPHA
    var_coefs    <- coefs %>% filter(term != "(Intercept)")
    all_sig      <- intercept_ok && length(var_coefs$p_value) > 0 &&
                     all(!is.na(var_coefs$p_value)) && all(var_coefs$p_value < SIG_ALPHA)

    formula_display <- paste(sapply(strsplit(gsub(".*~ (.*) \\+ \\(1.*", "\\1", fstr), " ")[[1]],
                                     function(w) if (w %in% c("+", "*")) w else display_name(w)),
                              collapse = " ")

    # Model-level columns (Formula/weight/AICc/DeltaAICc) are populated only
    # on the model's first row (the intercept) and left blank below it, so
    # the report reads as one block per model. The intercept is included as
    # its own "variable" row since its significance is checked anyway (see
    # all_sig below) and showing it makes each block self-contained.
    all_terms <- bind_rows(coefs[coefs$term == "(Intercept)", ], var_coefs)
    # N is a property of the phenotype's fitted data subset, not of any one
    # model (every candidate model for a phenotype is fit on the same
    # drop_na()-filtered rows) — repeated on every row rather than only the
    # first, since unlike Formula/weight/AICc it doesn't vary within a block.
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
      plotted[[as.character(rank)]] <- plot_response_surface(response, fstr, data, family, phenotype_label, out_path)
      # plot_response_lines() (the interaction-line-plot alternative) was
      # tried and dropped: with correlated predictors (common among
      # fog/BIO18/BIO15/CWD here), fixing one at an extreme percentile
      # while sweeping the other walks outside the observed data's convex
      # hull, so 2 of its 3 lines came out nearly or entirely masked away.
    }
  }
  invisible(plotted)
}

run_phenotype("Count_disk_flourets", "Count disk florets", "poisson",
              modern %>% drop_na(Count_disk_flourets))
run_phenotype("Count_ray_flourets", "Count ray florets", "poisson",
              modern %>% drop_na(Count_ray_flourets))
run_phenotype("Count_total_florets", "Count total florets", "poisson",
              modern %>% drop_na(Count_total_florets))
# Captured for the combined a/b/c/d figure at the very end of this script
# (panel_b = ray_prop, panel_a = Height_width) — see PHENOTYPE_DATA_PANELS.
PHENOTYPE_DATA_PANELS <- list()
ray_prop_plots <- run_phenotype("ray_prop", "Proportion ray florets", "beta",
              modern %>% drop_na(ray_prop))
if (length(ray_prop_plots) > 0) PHENOTYPE_DATA_PANELS$b <- ray_prop_plots[[1]]
run_phenotype("log(Height_in)", "Plant height (log inches)", "gaussian",
              modern %>% drop_na(Height_in))
run_phenotype("log(Width_in)", "Plant width (log inches)", "gaussian",
              modern %>% drop_na(Width_in))
height_width_plots <- run_phenotype("Height_width", "Height:width ratio", "gamma",
              modern %>% drop_na(Height_width))
if (length(height_width_plots) > 0) PHENOTYPE_DATA_PANELS$a <- height_width_plots[[1]]

# ============================================================
# 6. Wind turbine effect on ray floret count at Strauss
#
#    22 wind turbines were installed at Strauss (Tranquillon Mtn area)
#    around end of 2023. Naive per-locality before/after tests (as in the
#    original turbine_plot.R, now superseded/deleted) can't distinguish a
#    turbine-specific effect from a shared regional trend: the undergrad's
#    own stats_report.txt found Strauss's pre/post-2024 ray-count decline
#    (p=0.018-0.040 depending on the "before" window) was matched or
#    exceeded at VAFB_Arguello, a control site with no turbines (p=0.018
#    every window) — while Point_Conception, the other control, showed no
#    significant change (p=0.12-0.31). That pattern is exactly what a
#    formal period x treatment-group INTERACTION test is for: does the
#    before->after change at Strauss differ from the (pooled) controls'
#    change, rather than "is each site's own change significant".
#
#    Design decisions (2026-09-06, user-confirmed):
#      - Treatment = Locality (Strauss), not distance-to-turbine — wind
#        effects (fog/aridity) are directional (upwind/downwind), not
#        radially symmetric, and localities aren't a spatial continuum
#        (the nearest control site is genuinely far from Strauss), so a
#        distance gradient isn't the right model here.
#      - Ray floret count only, not disk — iNat photos (which supply most
#        of the pre-2024 "before" data) can't be used to count disk
#        florets, leaving only 2 non-NA disk-count specimens at Strauss in
#        the entire pre-2024 window; a period x treatment model for disk
#        count would be uninterpretable with that little data.
#      - Ray count fusion uses iNat WHEN AVAILABLE, falling back to
#        dissection — the reverse of this script's usual precedence
#        (dissection first elsewhere) — because a before/after design is
#        specifically vulnerable to the measurement method itself shifting
#        between periods; preferring one method uniformly where possible
#        keeps that consistent across the comparison. (The paired-
#        validation test in Section 1 found no significant difference
#        between the two methods either way, so this doesn't contradict
#        that result — it's a design choice, not a correction.)
#      - Env covariates use _germyr (that single specific year's own
#        conditions), not _modern (a 10-year rolling mean) — appropriate
#        here since the question is "did this specific year's conditions
#        or the turbines drive the change", not "how does the site's
#        typical climate relate to the trait". fog/BIO18/BIO15_germyr are
#        used (no CWD_germyr has been extracted).
#      - Before window: 2020-2023 vs. 2024 (the original script also
#        tried 2013-2023 and 2018-2023 as a robustness check; dropped
#        here per user decision — 2020-2023 is the one used).
# ============================================================

cat("\n=== Turbine effect: ray floret count, Strauss vs. controls ===\n\n")

turbine_localities <- c("Strauss", "VAFBArguello", "PointConception")

# iNat-preferred ray count fusion (opposite precedence from the rest of
# this script — see rationale above), built directly from pheno + inat
# rather than reusing `modern`, which uses dissection-preferred fusion.
turbine_base <- pheno %>%
  filter(!is.na(year), Locality %in% turbine_localities) %>%
  transmute(SampleID = recordNumber_CollectorNumber,
            year, Locality,
            Lat = `Lat._dd.dddd`, Long = Lon_dd.dddd,
            Count_ray_flourets = ifelse(!is.na(iNat_ray_count), iNat_ray_count, Count_ray_flourets)) %>%
  distinct()

turbine_inat_expansion <- inat %>%
  filter(!is.na(iNat_ray_count), !is.na(Year), !is.na(latitude), !is.na(longitude),
         normalize_locality(Locality) %in% turbine_localities,
         !(trimws(as.character(inat_Specimen)) %in% turbine_base$SampleID)) %>%
  transmute(SampleID = trimws(as.character(inat_Specimen)), year = Year,
            Locality = normalize_locality(Locality), Lat = latitude, Long = longitude,
            Count_ray_flourets = iNat_ray_count) %>%
  distinct()

turbine_base <- bind_rows(turbine_base, turbine_inat_expansion)
turbine_base$Locality <- factor(turbine_base$Locality, levels = turbine_localities)
turbine_base$treatment_group <- factor(ifelse(turbine_base$Locality == "Strauss", "treatment", "control"),
                                        levels = c("control", "treatment"))

# _germyr env covariates (that single specific collection year's own
# conditions, not the 10-year rolling _modern window used elsewhere)
env_germyr <- read.csv(env_file, stringsAsFactors = FALSE)
env_germyr$SampleID <- trimws(as.character(env_germyr$SampleID))
env_germyr <- distinct(env_germyr[, c("SampleID", "fog_germyr", "BIO18_germyr", "BIO15_germyr")])
turbine_base <- inner_join(turbine_base, env_germyr, by = "SampleID")
cat(sprintf("Turbine-analysis specimens (3 localities, any year, with germyr env data): %d\n",
            nrow(turbine_base)))

# --- descriptive plots (ported from turbine_plot.R; two bugs fixed: the
#     original used `data$Longitude`, a column that doesn't exist — the
#     real column is Lon_dd.dddd, so its map silently plotted 0 rows; and
#     an undefined variable `tr` in a plot title) ---

plot_period_boxplot <- function(data, cutoff_year, out_path) {
  d <- data %>% filter(year >= cutoff_year) %>%
    mutate(period = factor(ifelse(year < 2024, paste0(cutoff_year, "-2023"), "2024"),
                            levels = c(paste0(cutoff_year, "-2023"), "2024")))
  p <- ggplot(d, aes(x = period, y = Count_ray_flourets)) +
    geom_boxplot(fill = "gray80", color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.3, height = 0, size = 3, alpha = 0.4) +
    facet_wrap(~Locality, nrow = 1) +
    labs(x = NULL, y = "Ray floret count") +
    theme_bw(base_size = PLOT_BASE_SIZE) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = CONTOUR_LABEL_PT),
          axis.text  = element_text(size = CONTOUR_LABEL_PT * 0.7),
          strip.text = element_text(size = CONTOUR_LABEL_PT * 0.7))
  ggsave(out_path, p, width = 12, height = 6, dpi = 150)
}
plot_period_boxplot(turbine_base, 2020,
                     file.path(out_dir, "turbine_ray_before_after_2020.png"))

turbines <- read.csv("/private/groups/meyerlab/esteiger/pheno/Turbines.csv", stringsAsFactors = FALSE)
# Turbines were installed ~end of 2023, so only the 2024 facet should show
# them — giving `turbines` its own `period` column (rather than leaving it
# unfaceted, which would replicate the layer into every panel) restricts
# it to the panel where period=="2024".
turbines$period <- "2024"
map_data <- turbine_base %>% filter(Locality == "Strauss", year > 2020, !is.na(Long), !is.na(Lat)) %>%
  mutate(period = ifelse(year < 2024, "2021-2023", "2024"))
p_map <- ggplot(map_data, aes(x = Long, y = Lat)) +
  geom_jitter(width = 0.005, height = 0.005, aes(size = Count_ray_flourets), alpha = 0.3, color = "black") +
  geom_point(data = turbines, aes(x = Lon, y = Lat), shape = 13, color = "blue", size = 4, inherit.aes = FALSE) +
  facet_wrap(~period, nrow = 1) +
  scale_size(range = c(2, 7), name = "Ray count") +
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed() +
  theme_bw(base_size = PLOT_BASE_SIZE) +
  theme(axis.text = element_text(size = CONTOUR_LABEL_PT * 0.6, angle = 45),
        axis.title = element_text(size = CONTOUR_LABEL_PT * 0.8),
        strip.text = element_text(size = CONTOUR_LABEL_PT * 0.8),
        panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "turbine_map.png"), p_map, width = 9, height = 4.5, dpi = 150)

# --- formal difference-in-differences model, one per "before" window ---

turbine_report_rows <- list()
fit_turbine_window <- function(cutoff_year) {
  d <- turbine_base %>%
    filter(year >= cutoff_year) %>%
    mutate(period = factor(ifelse(year < 2024, "before", "2024"), levels = c("before", "2024"))) %>%
    drop_na(Count_ray_flourets, period, treatment_group, fog_germyr, BIO18_germyr, BIO15_germyr)
  d$Locality <- droplevels(d$Locality)
  d$fog_z   <- as.numeric(scale(d$fog_germyr))
  d$BIO18_z <- as.numeric(scale(d$BIO18_germyr))
  d$BIO15_z <- as.numeric(scale(d$BIO15_germyr))

  cat(sprintf("\n--- Before window: %d-2023 vs. 2024 (n=%d) ---\n", cutoff_year, nrow(d)))
  print(table(d$Locality, d$period))

  m <- fit_one("Count_ray_flourets ~ period * treatment_group + fog_z + BIO18_z + BIO15_z + (1 | Locality)",
               d, "poisson")
  if (is.null(m)) {
    warning("Turbine model failed to converge for window ", cutoff_year)
    return(invisible(NULL))
  }
  cc <- get_coefs(m)
  for (r in seq_len(nrow(cc))) {
    turbine_report_rows[[length(turbine_report_rows) + 1]] <<- data.frame(
      Window = sprintf("%d-2023 vs 2024", cutoff_year), N = nobs(m),
      Variable = cc$term[r], p_value = signif(cc$p_value[r], 2),
      Estimate_SE = sprintf("%.4f ± %.4f", cc$estimate[r], cc$std_error[r]),
      statistic = round(cc$statistic[r], 3)
    )
  }
  did_term <- "period2024:treatment_grouptreatment"
  cat(sprintf("  DiD estimate (period x treatment_group): est=%.4f ± %.4f, z=%.3f, p=%.3g%s\n",
              cc$estimate[cc$term == did_term], cc$std_error[cc$term == did_term],
              cc$statistic[cc$term == did_term], cc$p_value[cc$term == did_term],
              if (cc$p_value[cc$term == did_term] < SIG_ALPHA) "  *" else ""))
}
fit_turbine_window(2020)

turbine_report <- bind_rows(turbine_report_rows) %>%
  rename(`p-value` = p_value, `Estimate (β) ± std error` = Estimate_SE, `z-score` = statistic)
turbine_report_path <- file.path(out_dir, "turbine_report.csv")
write.csv(turbine_report, turbine_report_path, row.names = FALSE)
cat(sprintf("\nTurbine report written to %s\n", turbine_report_path))

# ============================================================
# 7. Mediation test: do Werner-fog/BIO18 trends over time explain each
#    floret phenotype's year-trend within the modern dataset?
#
#    Parallel to historic_data.R's mediation section (temporarily
#    disabled there while this version is built, 2026-09-06), but using
#    the MODIS-based Werner fog product (fog_modern — the germination-
#    year rolling window already used above, 2009-2024 coverage) instead
#    of Kling. Same three floret phenotypes as the historic version
#    (Count_disk_flourets, Count_ray_flourets, ray_prop); Height/Width/
#    ratio aren't tested here either, matching that script's scope.
#
#    Modern's collection years cluster heavily in 2020-2024 (rather than
#    spreading evenly over 2009-2024), so a secular trend has much less
#    time/spread to manifest than historic's ~65yr span — treat any
#    year-trend finding here as exploratory, same caveat as there.
# ============================================================

cat("\n=== Section 7: Mediation test (Werner fog) — do Fog/BIO18 trends over time explain each phenotype's year-trend? ===\n")

# 7a. Does each floret phenotype trend with collection year (locality-
#     adjusted), within the modern dataset?
cat("\n--- Year-trend models (trait ~ year, locality as random intercept) ---\n\n")
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
}
run_year_trend("Count_disk_flourets", "Count disk florets", "poisson", modern %>% drop_na(Count_disk_flourets))
run_year_trend("Count_ray_flourets", "Count ray florets", "poisson", modern %>% drop_na(Count_ray_flourets))
run_year_trend("Count_total_florets", "Count total florets", "poisson", modern %>% drop_na(Count_total_florets))
run_year_trend("ray_prop", "Proportion ray florets", "beta", modern %>% drop_na(ray_prop))

# 7b. Does Werner fog (or BIO18) itself trend over time? Uses
#     env_wi_range.csv's per-coordinate annual values (fog: 2000-2025;
#     BIO18: 1927-2024 — the full PRISM archive, same as historic_data.R
#     uses), locality-capped (15/coordinate, seed 42) per this project's
#     standard.
cat("\n--- Env-variable-over-time trend (env_wi_range.csv, locality-capped) ---\n\n")
range_env <- read.csv("/private/groups/meyerlab/esteiger/pheno/env_wi_range.csv", stringsAsFactors = FALSE)
set.seed(42)
range_capped <- do.call(rbind, lapply(split(range_env, range_env$Locality), function(d) {
  if (nrow(d) > 15) d[sample(nrow(d), 15), ] else d
}))
cat(sprintf("Locality-capped coordinates: %d (from %d total)\n", nrow(range_capped), nrow(range_env)))

pivot_env_long <- function(data, var_prefix) {
  cols <- grep(paste0("^", var_prefix, "_[0-9]{4}$"), names(data), value = TRUE)
  yrs  <- as.integer(sub(paste0("^", var_prefix, "_"), "", cols))
  data.frame(Locality = rep(data$Locality, length(cols)),
             Lat = rep(data$Lat, length(cols)), Long = rep(data$Long, length(cols)),
             year = rep(yrs, each = nrow(data)),
             value = unlist(data[, cols], use.names = FALSE))
}
run_env_trend <- function(var_prefix, var_label) {
  long <- pivot_env_long(range_capped, var_prefix) %>% drop_na(value)
  long$Locality <- droplevels(as.factor(long$Locality))
  long$year_z <- as.numeric(scale(long$year))
  m <- lmerTest::lmer(value ~ year_z + (1 | Locality), data = long)
  cc <- get_coefs(m) %>% filter(term == "year_z")
  cat(sprintf("  %s (n=%d obs, %d coords x %d yrs %d-%d): year beta = %.4f +/- %.4f, t = %.3f, p = %.3g%s\n",
              var_label, nobs(m), length(unique(paste(long$Lat, long$Long))), length(unique(long$year)),
              min(long$year), max(long$year), cc$estimate, cc$std_error, cc$statistic, cc$p_value,
              if (cc$p_value < SIG_ALPHA) "  *" else ""))
}
run_env_trend("fog", "Fog (Werner)")
run_env_trend("BIO18", "BIO18 (summer precip)")

# 7c. Mediation test proper: does adding Fog_z/BIO18_z (the Werner-based
#     _modern env variables already z-scored in Section 2) to the
#     year-trend model attenuate the year effect?
run_mediation <- function(response, label, family) {
  cat(sprintf("\n--- %s ---\n\n", label))
  med_data <- modern %>% drop_na(!!sym(response), year, Fog_z, BIO18_z)
  med_data$Locality <- droplevels(med_data$Locality)
  med_data$year_z <- as.numeric(scale(med_data$year))
  cat(sprintf("Full model (Fog + BIO18): n = %d specimens\n", nrow(med_data)))

  m_year_only <- fit_one(paste(response, "~ year_z + (1 | Locality)"), med_data, family)
  m_full      <- fit_one(paste(response, "~ year_z + Fog_z + BIO18_z + (1 | Locality)"), med_data, family)

  if (!is.null(m_year_only) && !is.null(m_full)) {
    cc_yo <- get_coefs(m_year_only) %>% filter(term == "year_z")
    cat(sprintf("  year alone (same n=%d, for fair comparison): beta=%.4f se=%.4f z=%.3f p=%.3g\n",
                nobs(m_year_only), cc_yo$estimate, cc_yo$std_error, cc_yo$statistic, cc_yo$p_value))
    cc_full <- get_coefs(m_full)
    for (r in seq_len(nrow(cc_full))) {
      cat(sprintf("  full model %-12s beta=%.4f se=%.4f z=%.3f p=%.3g%s\n",
                  paste0(cc_full$term[r], ":"), cc_full$estimate[r], cc_full$std_error[r],
                  cc_full$statistic[r], cc_full$p_value[r],
                  if (cc_full$p_value[r] < SIG_ALPHA) "  *" else ""))
    }
  } else {
    cat("  Full model (year + Fog + BIO18) failed to converge.\n")
  }
  invisible(NULL)
}
run_mediation("ray_prop", "Proportion ray florets", "beta")
run_mediation("Count_ray_flourets", "Count ray florets", "poisson")
run_mediation("Count_disk_flourets", "Count disk florets", "poisson")
run_mediation("Count_total_florets", "Count total florets", "poisson")

# ============================================================
# 8. Write report
# ============================================================

report <- bind_rows(report_rows) %>%
  rename(`Model weight` = Model_weight, `AICc` = AICc, `Delta AICc` = Delta_AICc,
         `p-value` = p_value, `Estimate (β) ± std error` = Estimate_SE, `z/t-score` = statistic)

report_path <- file.path(out_dir, "model_selection_report.csv")

# Multiple-comparisons caveat: with n_candidate_models candidate models fit
# per phenotype (all additive + pairwise-interaction combinations of the 4
# env variables) across 6 phenotypes, this report contains many p-values.
# None are corrected for multiple comparisons (across models or terms) —
# this is exploratory model selection by AICc, not a confirmatory,
# pre-registered set of hypothesis tests. Treat p < 0.05 here as
# hypothesis-generating, not as confirmed effects.
n_candidate_models <- length(build_formulas("y", z_vars))
caveat <- sprintf(
  "# Exploratory model selection: %d candidate models (all additive + pairwise-interaction combinations of %s) fit per phenotype, ranked by AICc. p-values are NOT corrected for multiple comparisons across models or terms — treat as hypothesis-generating, not confirmatory.",
  n_candidate_models, paste(names(env_vars), collapse = "/"))
writeLines(caveat, report_path)
suppressWarnings(write.table(report, report_path, sep = ",", row.names = FALSE, na = "",
                              append = TRUE, qmethod = "double"))
cat(sprintf("\nReport written to %s (%d rows, plus a multiple-comparisons caveat header)\n",
            report_path, nrow(report)))
