# phenotype_data.R
#
# Consolidated 2026-09-07 from modern_data.R + historic_data.R (kept as
# separate files no longer — this is the single source of truth for both
# analyses now). All modern-dataset analysis runs first, then all
# historic-dataset analysis, then a combined a/b/c/d summary figure at
# the very end (Section 11) built from panels captured along the way:
#   a: Height:width ratio response curve (modern)
#   b: Proportion ray florets response surface (modern)
#   c: BIO18-over-time pooled LOESS (historic)
#   d: Ray/disk floret counts over time (historic)
#
# Below this header: modern_data.R's original content, unchanged except
# for the run_phenotype() panel-capture additions (see Section 5) and
# the run_phenotype() same-N-across-candidates fix (see 2026-09-07
# comment there), then historic_data.R's original content, unchanged
# except for the Section 11 combined-figure addition at the very end.

# modern_data.R
#
# Reworked, self-contained version of the undergrad's modern-specimen GLM
# analyses (see gvtp_popgen/pheno/modern_models.R, heightmodels.R,
# updatedvarmodels.R for the original, ad hoc versions this replaces).
#
# For every phenotype record with a computable env window (effectively
# germination year 2009-2024 — see Section 2), fits all additive and
# interaction combinations of three env variables — fog and BIO18 (summer
# precip), each a rolling 20-year mean ending in (and including) the
# specimen's own germination year, not a fixed reference window (see
# extract_pheno_env_data.R Section 15); and CWD_25yr (climatic water
# deficit), a fixed-window (2000-2025) per-location mean, added
# 2026-10-22 after confirming it isn't excessively collinear with fog/
# BIO18 (Spearman rho 0.52 / -0.61, vs. -0.56 between fog and BIO18
# themselves) and costs no sample size (populated for every specimen
# already in the fog/BIO18-eligible set). BIO15 (precip seasonality) was
# dropped 2026-09-06 — collinear with fog (Spearman rho -0.84 to -0.89
# depending on window length) beyond this project's 0.7 threshold —
# against seven response variables:
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
# Wilcoxon test on double-measured specimens (moved to
# fire_turbine_analyses.R 2026-10-06, alongside the turbine analysis
# that also depends on iNat ray counts — see that script) showing no
# significant difference between the two methods.
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
  library(patchwork)  # side-by-side panel composition (Section 2c)
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
AXIS_LABELS <- c(fog_modern = "Percent hours with fog", BIO18_modern = "Summer precip. (mm)",
                  BIO15_modern = "Precip. seasonality", CWD_25yr = "Climatic water deficit (mm)")
PLOT_BASE_SIZE     <- 20
CONTOUR_LABEL_SIZE <- 7
# geom_label()'s `size` is in mm; convert to the pt size it actually
# renders at (via ggplot2's own mm-to-pt factor) so axis title/text can be
# set to match it exactly, rather than relying on theme_bw()'s default
# (smaller, relative-to-base-size) axis text scaling.
CONTOUR_LABEL_PT <- CONTOUR_LABEL_SIZE * ggplot2::.pt

# ============================================================
# 1. Ray-floret count: fuse dissection and iNat-photo counts
#
#     Fusing dissection-based Count_ray_flourets with iNat_ray_count
#     (for specimens with no physical collection) is justified by a
#     paired Wilcoxon signed-rank test on double-measured specimens —
#     moved to fire_turbine_analyses.R 2026-10-06 (run that script for
#     the validation stats; result: not significant, so the two are
#     treated as interchangeable here — dissection where available,
#     iNat count otherwise).
# ============================================================

pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
pheno$recordNumber_CollectorNumber <- trimws(as.character(pheno$recordNumber_CollectorNumber))
pheno$year <- suppressWarnings(as.numeric(pheno$year))
pheno$Count_ray_flourets <- suppressWarnings(as.numeric(pheno$Count_ray_flourets))
pheno$iNat_ray_count     <- suppressWarnings(as.numeric(pheno$iNat_ray_count))

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
env <- distinct(env[, c("SampleID", "fog_modern", "BIO18_modern", "BIO15_modern",
                        "fog_germyr", "BIO18_germyr", "CWD_25yr")])

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

# z-scored predictors used for model selection. BIO15 (BIO15_modern)
# dropped again 2026-09-07 after a re-test (2026-09-07, briefly re-added)
# confirmed the original 2026-09-06 finding: Fog x BIO15 collinearity
# (Spearman rho -0.84 at a 10yr rolling window, -0.89 at 20yr) makes any
# BIO15-containing model unstable — e.g. Height:width's top model became
# "Fog * BIO15" at weight=0.799, but its data points fall almost exactly
# along one diagonal line in Fog x BIO15 space (barely any real
# 2-D variation), so that "surface" was mostly extrapolation.
# CWD (CWD_25yr) added 2026-10-22: a fixed WY2000-2025 per-location mean
# (not germination-year-rolling like fog/BIO18 — no BCM data exists at
# the per-germination-year resolution fog/PRISM have), but populated for
# every specimen already in the fog/BIO18-eligible set, so it costs no
# sample size. Checked for the same instability BIO15 showed: CWD's
# correlation with fog (rho=0.52) and BIO18 (rho=-0.61) is no worse than
# fog-BIO18's own (-0.56), and in practice CWD never wins or meaningfully
# shifts the Fog+BIO18 coefficients for either focal phenotype (ray_prop,
# Height_width) — see modern_data_output/model_selection_report.csv.
env_vars <- c(Fog = "fog_modern", BIO18 = "BIO18_modern", CWD = "CWD_25yr")
for (nm in names(env_vars)) {
  modern[[paste0(nm, "_z")]] <- as.numeric(scale(modern[[env_vars[[nm]]]]))
}
z_vars <- paste0(names(env_vars), "_z")   # c("Fog_z","BIO18_z","CWD_z")

# ============================================================
# 2b. Locality (site) differences in modern floret counts and plant-size
#     metrics: Kruskal-Wallis omnibus test per phenotype, with Dunn's post
#     hoc (Bonferroni-corrected) pairwise comparisons for any phenotype
#     where the omnibus test is significant.
#
#     Ported and re-run here (2026-09-24) from the undergrad's original
#     stattests.R, which ran this on an older phenotype/env file and only
#     covered the floret-count phenotypes; extended here to also cover
#     Height_in/Width_in/Height_width, mirroring the GLMM section above.
#     Restricted to modern (year > 2000) specimens; CDFW and Dangermond
#     excluded (undersampled — n=20 and n=4 respectively — for a
#     fixed-effect locality comparison), exactly as in stattests.R.
# ============================================================

suppressMessages(library(FSA))   # dunnTest()

kw_dunn_data <- modern %>%
  filter(!Locality %in% c("CDFW", "Dangermond")) %>%
  mutate(Locality = droplevels(as.factor(Locality)))

cat(sprintf("\nSite-difference test dataset: %d records, %d localities (%s)\n",
            nrow(kw_dunn_data), nlevels(kw_dunn_data$Locality),
            paste(levels(kw_dunn_data$Locality), collapse = ", ")))

run_kw_dunn <- function(response, data, label) {
  d <- data %>% filter(!is.na(.data[[response]]))
  d$Locality <- droplevels(d$Locality)
  cat(sprintf("\n--- %s (n=%d, %d localities) ---\n", label, nrow(d), nlevels(d$Locality)))
  if (nlevels(d$Locality) < 2 || nrow(d) < 3) {
    cat("  Too few localities/observations, skipping.\n")
    return(NULL)
  }
  kw <- kruskal.test(reformulate("Locality", response), data = d)
  cat(sprintf("  Kruskal-Wallis: chi-sq = %.3f, df = %d, p = %.4g\n",
              kw$statistic, kw$parameter, kw$p.value))
  out <- list(kw = kw, dunn = NULL, n = nrow(d), label = label)
  if (kw$p.value < SIG_ALPHA) {
    dt <- dunnTest(d[[response]], d$Locality, method = "bonferroni")
    cat("  Significant -> Dunn's post hoc (Bonferroni):\n")
    print(dt$res)
    out$dunn <- dt$res
  }
  out
}

kw_dunn_results <- list(
  Count_ray_flourets  = run_kw_dunn("Count_ray_flourets",  kw_dunn_data, "Ray floret count"),
  Count_disk_flourets = run_kw_dunn("Count_disk_flourets", kw_dunn_data, "Disk floret count"),
  ray_prop            = run_kw_dunn("ray_prop",            kw_dunn_data, "Ray floret proportion"),
  Height_in           = run_kw_dunn("Height_in",           kw_dunn_data, "Plant height"),
  Width_in            = run_kw_dunn("Width_in",            kw_dunn_data, "Plant width"),
  Height_width        = run_kw_dunn("Height_width",        kw_dunn_data, "Height:width ratio")
)

kw_report_rows <- bind_rows(lapply(names(kw_dunn_results), function(nm) {
  r <- kw_dunn_results[[nm]]
  if (is.null(r)) return(NULL)
  data.frame(Phenotype = r$label, N = r$n, ChiSq = round(unname(r$kw$statistic), 3),
             df = unname(r$kw$parameter), p_value = signif(r$kw$p.value, 4))
}))
write.csv(kw_report_rows, file.path(out_dir, "locality_kw_report.csv"), row.names = FALSE)
cat(sprintf("\nWritten %s\n", file.path(out_dir, "locality_kw_report.csv")))

dunn_report_rows <- bind_rows(lapply(names(kw_dunn_results), function(nm) {
  r <- kw_dunn_results[[nm]]
  if (is.null(r) || is.null(r$dunn)) return(NULL)
  d <- r$dunn
  d$Phenotype <- r$label
  d
}))
if (nrow(dunn_report_rows) > 0) {
  write.csv(dunn_report_rows, file.path(out_dir, "locality_dunn_report.csv"), row.names = FALSE)
  cat(sprintf("Written %s\n", file.path(out_dir, "locality_dunn_report.csv")))
}

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
compute_contour_lines <- function(x_seq, y_seq, z, n_levels = 5, breaks = NULL) {
  if (all(is.na(z))) return(NULL)
  if (is.null(breaks)) breaks <- pretty(range(z, na.rm = TRUE), n = n_levels)
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

# contour_breaks: fixed contour levels (NULL = pretty()); label_size /
# point_size: contour-label (mm) and point sizes, so the combined figure
# (Section 12) can call this at its own smaller scale. The surface is shown
# in black/white (darker = higher predicted value; lightest shade kept well
# off white so it doesn't merge into the background) and is clipped to the
# observed range of the trait, as well as to the predictors' convex hull.
plot_response_surface <- function(response, formula_str, data, family, phenotype_label, out_path,
                                  contour_breaks = NULL, label_size = CONTOUR_LABEL_SIZE,
                                  point_size = 1.8) {
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

  # observed trait range on the response scale (gaussian models are fit on
  # log(response) but predictions are back-transformed)
  obs_trait <- eval(parse(text = response), data)
  if (family == "gaussian") obs_trait <- exp(obs_trait)
  trait_rng <- range(obs_trait, na.rm = TRUE)
  in_trait  <- function(pred) !is.na(pred) & pred >= trait_rng[1] & pred <= trait_rng[2]

  if (length(held_vars) == 0) {
    # Exactly 2 predictors — single panel, no held variables.
    grid <- base_grid
    grid$pred <- pred_fun(grid)
    grid$pred_masked <- ifelse(grid$in_hull & in_trait(grid$pred), grid$pred, NA)
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
    grid$pred_masked <- ifelse(grid$in_hull & in_trait(grid$pred), grid$pred, NA)
  }

  # Continuous raster + a few labeled contour lines, computed per facet
  # (each facet is its own prediction surface) — see compute_contour_lines().
  x_seq <- sort(unique(base_grid[[axis_vars[1]]]))
  y_seq <- sort(unique(base_grid[[axis_vars[2]]]))
  if (length(held_vars) == 0) {
    co <- compute_contour_lines(x_seq, y_seq, grid$pred_masked, breaks = contour_breaks)
    contour_lines  <- co$lines
    contour_labels <- co$labels
  } else {
    facet_groups <- levels(grid$facet_label)
    co_list <- lapply(facet_groups, function(fl) {
      sub <- grid[grid$facet_label == fl, ]
      co  <- compute_contour_lines(x_seq, y_seq, sub$pred_masked, breaks = contour_breaks)
      if (is.null(co)) return(NULL)
      co$lines$facet_label  <- factor(fl, levels = facet_groups)
      co$labels$facet_label <- factor(fl, levels = facet_groups)
      co
    })
    contour_lines  <- bind_rows(lapply(co_list, `[[`, "lines"))
    contour_labels <- bind_rows(lapply(co_list, `[[`, "labels"))
  }

  p <- ggplot() +
    geom_raster(data = grid, aes(x = .data[[axis_vars[1]]], y = .data[[axis_vars[2]]], fill = pred_masked)) +
    # Fog/BIO18/BIO15 are rolling-window means computed per germination
    # year, so many specimens sharing a GermYear (and similar coordinates)
    # land on the exact same point — a small jitter (~1.5% of each axis's
    # own range) breaks these apart so the true sample size is visible,
    # without meaningfully displacing where points fall relative to the
    # contour surface underneath them. Drawn before the contour lines/
    # labels below so those render on top and stay legible rather than
    # getting covered by a point.
    geom_point(data = data, aes(x = .data[[axis_vars[1]]], y = .data[[axis_vars[2]]]),
               shape = 21, fill = "white", color = "black", stroke = 0.4, size = point_size,
               position = position_jitter(width = diff(rng1) * 0.015, height = diff(rng2) * 0.015, seed = 42))
  if (!is.null(contour_lines) && nrow(contour_lines) > 0) {
    p <- p + geom_path(data = contour_lines, aes(x = x, y = y, group = id),
                        color = "white", linewidth = 0.4, alpha = 0.8)
  }
  if (!is.null(contour_labels) && nrow(contour_labels) > 0) {
    p <- p + geom_label(data = contour_labels, aes(x = x, y = y, label = signif(level, 3)),
                         size = label_size, linewidth = 0, fill = "white", alpha = 0.75,
                         color = "grey20", fontface = "bold")
  }
  p <- p +
    scale_fill_gradient(low = "grey78", high = "black", guide = "none", na.value = "transparent") +
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

RS_CALLS <- list()   # response-surface calls, by "response rank" (used by Section 12)
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
      RS_CALLS[[paste(response, rank)]] <<- list(formula = fstr, data = data, family = family,
                                                 label = phenotype_label)
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

# Panel a forced to the top-ranked "Fog + BIO18" model specifically (per
# user request 2026-09-07), bypassing the usual all-coefficients-significant
# plotting rule — that model's intercept isn't significant (p=0.24), which
# is why run_phenotype() above never plots it on its own (BIO18-alone,
# rank3, is what qualifies normally). Data subset matches this specific
# formula's own two predictors (Fog_z, BIO18_z) rather than the full
# (now 3-variable, CWD-included) z_vars, since this panel was never
# meant to additionally require CWD_25yr non-NA.
hw_panel_a_data <- modern %>% drop_na(Height_width, Fog_z, BIO18_z)
hw_panel_a_data$Locality <- droplevels(hw_panel_a_data$Locality)
PHENOTYPE_DATA_PANELS$a <- plot_response_surface(
  "Height_width", "Height_width ~ Fog_z + BIO18_z + (1 | Locality)",
  hw_panel_a_data, "gamma", "Height:width ratio",
  file.path(out_dir, "response_plot_Height_width_FogBIO18_forced.png"))

# ------------------------------------------------------------
# Panels B and C of the 6-panel figure (Section 12), drawn here while the
# modern section's plot_response_surface()/globals are still in effect
# (the historic section below redefines both). Same models as above, at
# the combined figure's scale: FIG_PT text, smaller labels/points.
#   B: Height:width ratio (Fog + BIO18), contour lines at 1-5
#   C: ray floret proportion (top significant model)
# ------------------------------------------------------------
FIG_PT       <- 7                  # all text in the 6-panel figure (pt)
FIG_THEME    <- theme(axis.title = element_text(size = FIG_PT),
                      axis.text  = element_text(size = FIG_PT - 0.5),
                      axis.ticks = element_line(linewidth = 0.3),
                      axis.ticks.length = unit(2, "pt"),
                      panel.border = element_rect(fill = NA, color = "grey20", linewidth = 0.4),
                      plot.margin = margin(1, 2, 1, 1))
PHENOTYPE_FIG_PANELS <- list()
PHENOTYPE_FIG_PANELS$B <- plot_response_surface(
  "Height_width", "Height_width ~ Fog_z + BIO18_z + (1 | Locality)",
  hw_panel_a_data, "gamma", "Height:width ratio",
  file.path(out_dir, "fig_panel_B_height_width.png"),
  contour_breaks = 1:5, label_size = (FIG_PT - 0.5) / ggplot2::.pt, point_size = 0.8) + FIG_THEME
rp_key  <- paste("ray_prop", names(ray_prop_plots)[1])
rp_call <- RS_CALLS[[rp_key]]
PHENOTYPE_FIG_PANELS$C <- plot_response_surface(
  "ray_prop", rp_call$formula, rp_call$data, rp_call$family, rp_call$label,
  file.path(out_dir, "fig_panel_C_ray_prop.png"),
  label_size = (FIG_PT - 0.5) / ggplot2::.pt, point_size = 0.8) + FIG_THEME
PHENOTYPE_FIG_PANELS$C_data <- rp_call$data
cat(sprintf("Figure panels B/C drawn (ray_prop model: %s)\n", rp_call$formula))

# ============================================================
# 5b. Height:width ratio and ray_prop, using germination-year (single
#     collection-year) env values instead of the 20yr rolling-window
#     means used above — same two variables (Fog, BIO18), same two
#     phenotypes, different source columns (fog_germyr/BIO18_germyr vs.
#     fog_modern/BIO18_modern). Only the top (rank1) model per phenotype
#     is checked, per user request (2026-09-11) — printed to console
#     regardless of outcome; plotted and added to the report only if the
#     intercept AND every coefficient are significant, same rule as
#     everywhere else in this script.
# ============================================================

env_vars_germyr <- c(Fog = "fog_germyr", BIO18 = "BIO18_germyr")
for (nm in names(env_vars_germyr)) {
  modern[[paste0(nm, "_germyr_z")]] <- as.numeric(scale(modern[[env_vars_germyr[[nm]]]]))
}
z_vars_germyr <- paste0(names(env_vars_germyr), "_germyr_z")   # c("Fog_germyr_z","BIO18_germyr_z")

run_phenotype_germyr <- function(response, phenotype_label, family, data) {
  cat(sprintf("\n=== %s, germyr env vars (%s) ===\n", phenotype_label, response))
  data <- data %>% drop_na(all_of(z_vars_germyr))
  data$Locality <- droplevels(data$Locality)
  fit <- fit_all(response, data, family, vars = z_vars_germyr)
  if (length(fit$models) == 0) {
    warning("No candidate models converged for ", response, " (germyr) — skipping.")
    return(invisible(NULL))
  }
  ranked <- aicc_rank(fit$models) %>% slice_head(n = 1)
  idx    <- ranked$index[1]
  model  <- fit$models[[idx]]
  fstr   <- fit$formulas[[idx]]
  coefs  <- get_coefs(model)

  intercept_ok <- coefs$p_value[coefs$term == "(Intercept)"] < SIG_ALPHA
  var_coefs    <- coefs %>% filter(term != "(Intercept)")
  all_sig      <- intercept_ok && length(var_coefs$p_value) > 0 &&
                   all(!is.na(var_coefs$p_value)) && all(var_coefs$p_value < SIG_ALPHA)

  cat(sprintf("  Top model: %s  (n=%d, weight=%.3f, AICc=%.2f)\n",
              fstr, nobs(model), ranked$Weight[1], ranked$AICc[1]))
  for (r in seq_len(nrow(coefs))) {
    cat(sprintf("    %-12s beta=%.4f se=%.4f z=%.3f p=%.3g%s\n",
                paste0(coefs$term[r], ":"), coefs$estimate[r], coefs$std_error[r],
                coefs$statistic[r], coefs$p_value[r],
                if (coefs$p_value[r] < SIG_ALPHA) "  *" else ""))
  }
  cat(sprintf("  Intercept + all coefficients significant: %s\n", if (all_sig) "YES" else "no"))

  if (!all_sig) return(invisible(NULL))

  cat("  -> plotting and adding to model_selection_report.csv.\n")
  all_terms <- bind_rows(coefs[coefs$term == "(Intercept)", ], var_coefs)
  n_obs <- nobs(model)
  disp_map <- setNames(names(env_vars_germyr), z_vars_germyr)
  formula_display <- paste(sapply(strsplit(gsub(".*~ (.*) \\+ \\(1.*", "\\1", fstr), " ")[[1]],
                                   function(w) if (w %in% c("+", "*")) w else str_replace_all(w, disp_map)),
                            collapse = " ")
  for (r in seq_len(nrow(all_terms))) {
    term_display <- if (all_terms$term[r] == "(Intercept)") "(Intercept)" else
      str_replace_all(all_terms$term[r], disp_map) %>% str_replace_all(":", " x ")
    report_rows[[length(report_rows) + 1]] <<- data.frame(
      Phenotype    = paste(phenotype_label, "(germyr env)"),
      N            = n_obs,
      Formula      = if (r == 1) formula_display else "",
      Model_weight = if (r == 1) round(ranked$Weight[1], 3) else NA,
      AICc         = if (r == 1) round(ranked$AICc[1], 2) else NA,
      Delta_AICc   = if (r == 1) 0 else NA,
      Variable     = term_display,
      p_value      = signif(all_terms$p_value[r], 2),
      Estimate_SE  = sprintf("%.3f ± %.3f", all_terms$estimate[r], all_terms$std_error[r]),
      statistic    = round(all_terms$statistic[r], 3)
    )
  }

  # plot_response_surface() reads env_vars/z_vars/AXIS_LABELS as globals —
  # swap them in temporarily so it can be reused as-is for this germyr
  # variable set, then restore, so the main (rolling-window) analysis
  # above is unaffected by this section.
  old_env_vars <- env_vars; old_z_vars <- z_vars; old_axis_labels <- AXIS_LABELS
  env_vars <<- env_vars_germyr
  z_vars   <<- z_vars_germyr
  AXIS_LABELS <<- c(AXIS_LABELS, fog_germyr = "Fog", BIO18_germyr = "Warmest quarter precip.")
  out_path <- file.path(out_dir, sprintf("response_plot_%s_germyr.png", response))
  p <- plot_response_surface(response, fstr, data, family, phenotype_label, out_path)
  env_vars <<- old_env_vars; z_vars <<- old_z_vars; AXIS_LABELS <<- old_axis_labels
  invisible(p)
}

run_phenotype_germyr("Height_width", "Height:width ratio", "gamma",
                      modern %>% drop_na(Height_width))
run_phenotype_germyr("ray_prop", "Proportion ray florets", "beta",
                      modern %>% drop_na(ray_prop))

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

# ============================================================
# ======================  HISTORIC DATASET  ====================
# ============================================================

# historic_data.R
#
# Companion to modern_data.R, following the undergrad's histenvvarsmodels.R
# in spirit (long-window fog via Kling instead of MODIS/Werner). Uses
# three _historic environmental variables — BIO18 (summer precip), Kling
# fog (the real observed 1950-2014 fog series; the analog of "fog" here,
# since MODIS FLCC only starts 2000), and CWD_historic (climatic water
# deficit, added 2026-10-26) — each a rolling 20-year mean ending in (and
# including) the specimen's own germination year, not a fixed reference
# window (see extract_pheno_env_data.R Section 15). BIO15 (precip
# seasonality) is excluded here for the same reason as the modern
# section above (collinear with fog beyond this project's 0.7
# threshold) — applied to the historic candidate set 2026-10-27 per
# user instruction, not independently re-tested here. Kling's 1950-2014
# coverage caps
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
# test as the modern section above (moved to fire_turbine_analyses.R
# 2026-10-06 — never year-restricted, so one run covers both sections).
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
  library(Kendall)    # non-parametric corroboration of run_env_trend()'s LMER trends
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
# 1. Ray-floret count: fuse dissection and iNat-photo counts
#    (identical logic/result to the modern section above — see there,
#    and fire_turbine_analyses.R for the validating Wilcoxon test)
# ============================================================

pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
pheno$recordNumber_CollectorNumber <- trimws(as.character(pheno$recordNumber_CollectorNumber))
pheno$year <- suppressWarnings(as.numeric(pheno$year))
pheno$Count_ray_flourets <- suppressWarnings(as.numeric(pheno$Count_ray_flourets))
pheno$iNat_ray_count     <- suppressWarnings(as.numeric(pheno$iNat_ray_count))

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
env <- distinct(env[, c("SampleID", "BIO18_historic", "BIO15_historic", "Kling_historic", "CWD_historic")])

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
historic$Source <- "dissection"

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
inat_expansion$Source <- "iNat"
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

# --- Does ray floret count differ systematically between dissection- and
# iNat-photo-derived counts? The ray/disk-over-time figure (Section 3
# below) pools both sources into one "Ray" series, and iNat records are
# concentrated in recent years (citizen-science platform), so a counting-
# method artifact here could masquerade as (or mask) a real temporal trend
# in that figure. Tests Count_ray_flourets ~ Source, locality-adjusted,
# same poisson-with-overdispersion-refit convention as run_year_trend()
# below (fit_one() auto-refits to negative binomial if overdispersed).
cat("\n--- Ray floret count: dissection vs. iNat-photo counting method ---\n\n")
source_check_data <- historic %>% drop_na(Count_ray_flourets, Source)
source_check_data$Locality <- droplevels(as.factor(source_check_data$Locality))
m_source <- fit_one("Count_ray_flourets ~ Source + (1 | Locality)", source_check_data, "poisson")
if (is.null(m_source)) {
  warning("Dissection-vs-iNat ray count model failed to converge")
} else {
  cc_source <- get_coefs(m_source) %>% filter(grepl("^Source", term))
  cat(sprintf("  Ray count, iNat vs. dissection (n=%d, %d dissection / %d iNat): beta = %.4f +/- %.4f, z/t = %.3f, p = %.3g%s\n",
              nobs(m_source), sum(source_check_data$Source == "dissection"),
              sum(source_check_data$Source == "iNat"),
              cc_source$estimate, cc_source$std_error, cc_source$statistic, cc_source$p_value,
              if (cc_source$p_value < SIG_ALPHA) "  *  -- counting methods differ, see caveat above" else "  (no evidence of a counting-method difference)"))
  write.csv(data.frame(N_obs = nobs(m_source),
                       N_dissection = sum(source_check_data$Source == "dissection"),
                       N_iNat = sum(source_check_data$Source == "iNat"),
                       `Estimate (beta) +/- std error` = sprintf("%.4f +/- %.4f", cc_source$estimate, cc_source$std_error),
                       `z/t-score` = round(cc_source$statistic, 3), `p-value` = signif(cc_source$p_value, 2),
                       check.names = FALSE),
            file.path(out_dir, "ray_count_dissection_vs_inat_report.csv"), row.names = FALSE)
}

# z-scored predictors used for model selection. BIO15 (BIO15_historic)
# excluded here for the same reason as the modern section above
# (2026-09-06 finding: Fog x BIO15 collinearity makes any BIO15-
# containing model unstable) -- not re-tested per-phenotype here, just
# applying the modern section's established exclusion to the historic
# candidate set as well, per user instruction 2026-10-27.
env_vars <- c(Fog = "Kling_historic", BIO18 = "BIO18_historic", CWD = "CWD_historic")
for (nm in names(env_vars)) {
  historic[[paste0(nm, "_z")]] <- as.numeric(scale(historic[[env_vars[[nm]]]]))
}
z_vars <- paste0(names(env_vars), "_z")   # c("Fog_z","BIO18_z","CWD_z")

# ============================================================
# 2b. Distributions of raw phenotype measurements (all specimens,
#     modern + historic + iNat-only), faceted by locality
#
#     Sanity-check plot: disk/ray floret count and ray proportion are
#     built directly from the raw `pheno`/`inat` tables (not `modern`/
#     `historic`, which are each restricted by an inner_join against a
#     rolling-window env-coverage window — irrelevant to whether a
#     floret count itself is trustworthy, and previously found to
#     silently drop 15 of 21 iNat-only specimens from this plot for
#     that reason, 2026-10-07). This recovers every specimen with a
#     disk or ray count, full 1927-2024 archive, no env-coverage gate.
#     Plant height/width/ratio, by contrast, were only ever recorded on
#     the 2024 modern collection (confirmed directly from the raw file
#     — all 115 non-NA Height_in rows are year==2024, and `modern`'s own
#     env-coverage window doesn't exclude any of them — not an
#     env-coverage artifact), so those three panels still show
#     `modern`'s specimens. Locality names relabeled to match this
#     project's other plotting scripts (plot_env_time.R's
#     cluster_labels); one row per locality (facet_wrap ncol=1), with
#     the two row-blocks' heights scaled to the number of localities
#     each actually has data for, so the shorter plant-size block isn't
#     stretched to match the taller floret-count block.
# ============================================================

LOCALITY_LABELS <- c(VAFBLionHead = "Lion's Head", Strauss = "Tranquillon Mtn",
                      VAFBArguello = "Pt. Arguello", Dangermond = "Dangermond",
                      PointConception = "Pt. Conception", Agujas = "Agujas",
                      CDFW = "Gaviota")
relabel_locality <- function(x) {
  x <- as.character(x)
  lab <- unname(LOCALITY_LABELS[x])
  out <- ifelse(is.na(lab), x, lab)
  lvls <- unique(c(unname(LOCALITY_LABELS), sort(setdiff(out, unname(LOCALITY_LABELS)))))
  factor(out, levels = lvls)
}

raw_floret <- pheno %>%
  transmute(SampleID = recordNumber_CollectorNumber, Locality = normalize_locality(Locality),
            Count_disk_flourets = suppressWarnings(as.numeric(Count_disk_flourets)),
            Count_ray_flourets) %>%
  distinct()

# Requires real coordinates (!is.na(latitude/longitude)) — matches modern's/
# historic's own iNat-expansion filter; without it, ~15 iNat rows with a
# ray count but no locality AND no coordinates at all (unplaceable, not
# just missing a locality label) were getting pulled in as a spurious
# "NA" locality facet (caught 2026-10-07 on first attempt at this fix).
raw_floret_inat <- inat %>%
  filter(!is.na(iNat_ray_count), !is.na(latitude), !is.na(longitude),
         !(SampleID %in% raw_floret$SampleID)) %>%
  transmute(SampleID, Locality = normalize_locality(Locality),
            Count_disk_flourets = NA_real_, Count_ray_flourets = iNat_ray_count) %>%
  distinct()

floret_combined <- bind_rows(raw_floret, raw_floret_inat)

# Backfill a still-missing Locality from inatrays.csv when available (same
# step modern/historic each do) — e.g. specimen 233000618, blank Locality
# in pheno_file, backfilled to "Coastal_Lompoc" from its iNat record.
na_loc <- is.na(floret_combined$Locality)
if (any(na_loc)) {
  inat_loc_by_id <- setNames(normalize_locality(inat$Locality), inat$SampleID)
  floret_combined$Locality[na_loc] <- unname(inat_loc_by_id[floret_combined$SampleID[na_loc]])
}

floret_combined <- floret_combined %>%
  distinct(SampleID, .keep_all = TRUE) %>%
  mutate(ray_prop = Count_ray_flourets / (Count_disk_flourets + Count_ray_flourets),
         Locality = relabel_locality(Locality))
n_inat_only <- nrow(raw_floret_inat)
cat(sprintf("Combined floret-count dataset for distribution plot: %d specimens (%d from physical collections, %d iNat-only)\n",
            nrow(floret_combined), nrow(floret_combined) - n_inat_only, n_inat_only))

modern_relabel <- modern %>% mutate(Locality = relabel_locality(Locality))

dist_hist <- function(data, var, label) {
  d <- data %>% transmute(value = .data[[var]], Locality) %>% drop_na(value)
  ggplot(d, aes(x = value)) +
    geom_histogram(fill = "black", color = "black") +
    facet_wrap(~Locality, ncol = 1, scales = "free_y") +
    labs(x = label, y = "Count") +
    theme_bw(base_size = 13) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
          strip.text = element_text(size = 9))
}

p_dist_disk   <- dist_hist(floret_combined, "Count_disk_flourets", "Disk floret count")
p_dist_ray    <- dist_hist(floret_combined, "Count_ray_flourets", "Ray floret count")
p_dist_rayp   <- dist_hist(floret_combined, "ray_prop", "Ray proportion")
p_dist_height <- dist_hist(modern_relabel, "Height_in", "Plant height (in)")
p_dist_width  <- dist_hist(modern_relabel, "Width_in", "Plant width (in)")
p_dist_hw     <- dist_hist(modern_relabel, "Height_width", "Height:width ratio")

n_locs_present <- function(data, var) n_distinct(data$Locality[!is.na(data[[var]])])
n_rows_top    <- max(sapply(c("Count_disk_flourets", "Count_ray_flourets", "ray_prop"),
                             n_locs_present, data = floret_combined))
n_rows_bottom <- max(sapply(c("Height_in", "Width_in", "Height_width"),
                             n_locs_present, data = modern_relabel))

p_dist_combined <- (p_dist_disk | p_dist_ray | p_dist_rayp) / (p_dist_height | p_dist_width | p_dist_hw) +
  plot_layout(heights = c(n_rows_top, n_rows_bottom))

PER_ROW_IN <- 1.3
ggsave(file.path(out_dir, "phenotype_distributions.png"), p_dist_combined,
       width = 10, height = PER_ROW_IN * (n_rows_top + n_rows_bottom) + 1, dpi = 150)
cat(sprintf("Phenotype distribution plot (faceted by locality) written to %s\n",
            file.path(out_dir, "phenotype_distributions.png")))

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

# Kling fog frequency, same pooled-LOESS treatment, added to panel c on
# a secondary (right-hand) axis — its scale (~0.05-0.6) is entirely
# different from BIO18's (~1-4.5), so its fitted values are linearly
# rescaled onto BIO18's axis for plotting, with sec_axis() applying the
# inverse transform back to Kling's own units for the right-hand labels.
#
# Restored 2026-10-07 (per user request) after being swapped out for
# BIO15 on 2026-09-24 — back to Kling on the right axis, BIO18 (summer
# precip) on the left, as in the original bio18_kling_over_time.png.
kling_long <- env_range %>%
  select(Locality, Lat, Long, matches("^kling_fog_[0-9]{4}$")) %>%
  pivot_longer(matches("^kling_fog_[0-9]{4}$"), names_to = "year", values_to = "value") %>%
  mutate(year = as.integer(sub("kling_fog_", "", year))) %>%
  drop_na(value) %>%
  # Kling's raw value is a 0-1 fraction; rescaled to a true percentage here
  # (display only — the underlying GLMM/trend-test variable, Kling_historic/
  # kling_fog, is untouched) so the axis and its "% fog" label agree, matching
  # fog_modern (Werner FLCC), which is already 0-100.
  mutate(value = value * 100)
# span=0.2 matches plot_env_time.R's own var_spans[["kling_fog"]] choice
# there (less smoothing than BIO18's 0.5 — Kling is a smoother, more
# autocorrelated physical process than storm-driven precipitation, so it
# tolerates less smoothing without getting noisy; BIO18 needs the extra
# smoothing specifically because it's zero-inflated with sparse spikes).
kling_loess <- trim_loess_bio18(kling_long$year, kling_long$value, span = 0.2)

FOG_LINE_COLOR <- "grey70"
FOG_RIBBON_COLOR <- "grey80"
kling_range <- range(kling_loess$fitted)
bio18_range <- range(c(bio18_loess$ymin, bio18_loess$ymax))
kb_slope <- diff(bio18_range) / diff(kling_range)
kb_intercept <- bio18_range[1] - kb_slope * kling_range[1]
kling_loess$fitted_scaled <- kb_intercept + kb_slope * kling_loess$fitted
kling_loess$ymin_scaled   <- kb_intercept + kb_slope * kling_loess$ymin
kling_loess$ymax_scaled   <- kb_intercept + kb_slope * kling_loess$ymax

p_bio18 <- ggplot(bio18_loess, aes(x = year, y = fitted)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4) +
  geom_ribbon(data = kling_loess, aes(x = year, y = fitted_scaled, ymin = ymin_scaled, ymax = ymax_scaled),
              fill = FOG_RIBBON_COLOR, alpha = 0.4) +
  geom_line(color = "black", linewidth = 0.6) +
  geom_line(data = kling_loess, aes(x = year, y = fitted_scaled), color = FOG_LINE_COLOR, linewidth = 1) +
  scale_x_continuous(breaks = c(1950, 1970, 1990, 2010)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14)),
                      sec.axis = sec_axis(transform = ~ (. - kb_intercept) / kb_slope,
                                          name = "Percent hours with fog")) +
  labs(x = "Year", y = "Summer precip. (mm)") +
  theme_bw(base_size = PLOT_BASE_SIZE) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = RD_TIME_FONT_PT),
        axis.text  = element_text(size = RD_TIME_FONT_PT * 0.75),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.title.y.right = element_text(size = RD_TIME_FONT_PT, color = FOG_LINE_COLOR),
        axis.text.y.right  = element_text(size = RD_TIME_FONT_PT * 0.75, color = FOG_LINE_COLOR),
        plot.margin = margin(t = 5, r = 8, b = 2, l = 5))

ggsave(file.path(out_dir, "bio18_kling_over_time.png"), p_bio18, width = 6.5, height = 6.5, dpi = 150)
cat(sprintf("Panel c (BIO18/Kling fog over time), standalone, written to %s\n",
            file.path(out_dir, "bio18_kling_over_time.png")))

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
  scale_x_continuous(breaks = c(1950, 1970, 1990, 2010)) +
  labs(x = "Year", y = "Floret count", shape = NULL, color = NULL) +
  theme_bw(base_size = PLOT_BASE_SIZE) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = RD_TIME_FONT_PT),
        axis.text  = element_text(size = RD_TIME_FONT_PT * 0.75),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 4)),
        legend.text = element_text(size = RD_TIME_FONT_PT * 0.75),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
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
    geom_raster(data = grid, aes(x = .data[[axis_vars[1]]], y = .data[[axis_vars[2]]], fill = pred_masked)) +
    # Fog/BIO18/BIO15 are rolling-window means computed per germination
    # year, so many specimens sharing a GermYear (and similar coordinates)
    # land on the exact same point — a small jitter (~1.5% of each axis's
    # own range) breaks these apart so the true sample size is visible,
    # without meaningfully displacing where points fall relative to the
    # contour surface underneath them. Drawn before the contour lines/
    # labels below so those render on top and stay legible rather than
    # getting covered by a point.
    geom_point(data = data, aes(x = .data[[axis_vars[1]]], y = .data[[axis_vars[2]]]),
               shape = 21, fill = "white", color = "black", stroke = 0.4, size = 1.8,
               position = position_jitter(width = diff(rng1) * 0.015, height = diff(rng2) * 0.015, seed = 42))
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
# Fitted model + the year mean/sd used to build year_z, captured per
# response so panel d (below) can plot this GLMM's own population-level
# prediction instead of a naive unadjusted lm() smooth.
year_trend_models <- list()
run_year_trend <- function(response, phenotype_label, family, data) {
  data <- data %>% drop_na(year)
  data$Locality <- droplevels(data$Locality)
  year_mean <- mean(data$year); year_sd <- sd(data$year)
  data$year_z <- (data$year - year_mean) / year_sd
  m <- fit_one(paste(response, "~ year_z + (1 | Locality)"), data, family)
  if (is.null(m)) {
    warning("Year-trend model failed to converge for ", response)
    return(invisible(NULL))
  }
  year_trend_models[[response]] <<- list(model = m, year_mean = year_mean, year_sd = year_sd,
                                          year_range = range(data$year))
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

# --- Rebuild panel d (combined-figure version only — the standalone
#     ray_disk_over_time.png above keeps its original naive lm smooth)
#     using each phenotype's own year-trend GLMM (fit just above) for
#     the trend line/CI, rather than an unadjusted lm() that ignores
#     the (1|Locality) structure the actual reported p-values come
#     from. Population-level (re.form=NA) prediction on the link scale,
#     back-transformed; SE from the fixed-effects vcov only (ignoring
#     between-locality variance, standard for this kind of "average
#     locality" band) since glmer's predict() has no built-in se.fit.
predict_year_trend_glmm <- function(fit) {
  m <- fit$model
  year_seq <- seq(fit$year_range[1], fit$year_range[2], length.out = 100)
  newdata <- data.frame(year_z = (year_seq - fit$year_mean) / fit$year_sd)
  if (inherits(m, "glmmTMB")) {
    pred <- predict(m, newdata = newdata, type = "link", se.fit = TRUE)
    link_fit <- pred$fit; link_se <- pred$se.fit
  } else {
    X <- model.matrix(~year_z, newdata)
    beta <- lme4::fixef(m)
    link_fit <- as.vector(X %*% beta)
    V <- as.matrix(vcov(m))
    link_se <- sqrt(rowSums((X %*% V) * X))
  }
  data.frame(year = year_seq, fit = exp(link_fit),
             ymin = exp(link_fit - 1.96 * link_se), ymax = exp(link_fit + 1.96 * link_se))
}

if (all(c("Count_disk_flourets", "Count_ray_flourets") %in% names(year_trend_models))) {
  disk_pred <- predict_year_trend_glmm(year_trend_models[["Count_disk_flourets"]])
  disk_pred$metric <- "Disk"
  ray_pred <- predict_year_trend_glmm(year_trend_models[["Count_ray_flourets"]])
  ray_pred$metric <- "Ray"
  time_pred <- bind_rows(disk_pred, ray_pred) %>%
    mutate(metric = factor(metric, levels = c("Disk", "Ray")))

  p_time_glmm <- ggplot(time_data, aes(x = year, y = value, shape = metric, color = metric)) +
    geom_point(position = position_jitter(width = 0.3, height = 0.3, seed = 42),
               alpha = 0.85, size = 2.3) +
    geom_ribbon(data = time_pred, aes(x = year, y = fit, ymin = ymin, ymax = ymax, group = metric),
                inherit.aes = FALSE, fill = "grey70", alpha = 0.4) +
    geom_line(data = time_pred, aes(x = year, y = fit, group = metric), inherit.aes = FALSE,
              color = "black", linewidth = 0.6) +
    scale_shape_manual(values = c("Disk" = 16, "Ray" = 17)) +
    scale_color_manual(values = c("Disk" = "grey65", "Ray" = "black")) +
    scale_x_continuous(breaks = c(1950, 1970, 1990, 2010)) +
    labs(x = "Year", y = "Floret count", shape = NULL, color = NULL) +
    theme_bw(base_size = PLOT_BASE_SIZE) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = RD_TIME_FONT_PT),
          axis.text  = element_text(size = RD_TIME_FONT_PT * 0.75),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_text(margin = margin(t = 4)),
          legend.text = element_text(size = RD_TIME_FONT_PT * 0.75),
          legend.title = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
          plot.margin = margin(t = 5, r = 8, b = 2, l = 5))
  PHENOTYPE_DATA_PANELS$d <- p_time_glmm
  cat("Panel d rebuilt using the year-trend GLMMs' own population-level predictions.\n")
} else {
  warning("Disk/ray year-trend models not both available — panel d keeps the naive lm smooth.")
}

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

# --- Non-parametric corroboration: Mann-Kendall test on each variable's
# locality-capped annual mean series (same range_capped set run_env_trend()
# uses above, same 15/coordinate cap, same seed -- so this test and the
# LMER trend above are drawn from identical input data, just a different
# statistic). Each variable's own full available year range is used (no
# year-range restriction is applied here, or anywhere above -- e.g. Kling
# fog's annual layers only exist for 1950-2014 in the first place, so that
# is simply its full record, not a chosen subset; BIO18/BIO15 use their
# full 1927-2024 PRISM archive).
mk_rows <- list()
run_mann_kendall <- function(var_prefix, var_label) {
  long <- pivot_env_long(range_capped, var_prefix) %>% drop_na(value)
  annual_mean <- long %>% group_by(year) %>% summarise(value = mean(value), .groups = "drop") %>%
    arrange(year)
  mk <- MannKendall(annual_mean$value)
  cat(sprintf("  %s (annual means, %d yrs %d-%d, %d coords): Mann-Kendall tau = %.4f, p = %.3g%s\n",
              var_label, nrow(annual_mean), min(annual_mean$year), max(annual_mean$year),
              length(unique(paste(long$Lat, long$Long))), mk$tau, mk$sl,
              if (mk$sl < SIG_ALPHA) "  *" else ""))
  mk_rows[[length(mk_rows) + 1]] <<- data.frame(
    Variable = var_label, N_years = nrow(annual_mean),
    Year_range = sprintf("%d-%d", min(annual_mean$year), max(annual_mean$year)),
    Tau = round(mk$tau, 4), `p-value` = signif(mk$sl, 2), check.names = FALSE)
}
cat("\n--- Mann-Kendall corroboration (annual means, locality-capped) ---\n\n")
run_mann_kendall("kling_fog", "Fog (Kling)")
run_mann_kendall("BIO18", "BIO18 (summer precip)")
run_mann_kendall("BIO15", "BIO15 (precip seasonality)")

mann_kendall_report <- bind_rows(mk_rows)
mann_kendall_path <- file.path(out_dir, "mann_kendall_report.csv")
write.csv(mann_kendall_report, mann_kendall_path, row.names = FALSE)
cat(sprintf("\nMann-Kendall report written to %s\n", mann_kendall_path))

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

  # --- Phenotype-extreme silhouettes (transparent PNGs supplied by user,
  #     2026-09-11) placed in the corners of panels a/b. annotation_custom()
  #     positions in DATA units, so each image's on-plot footprint is sized
  #     from the panel's own data range (height = a fraction of the y-range)
  #     with width back-solved from the image's true pixel aspect ratio and
  #     an estimated panel physical aspect ratio (panel render ~4.2in wide x
  #     6.0in tall, from the 18x7in/4-panel layout minus axis margins) —
  #     approximate, since patchwork doesn't expose exact per-panel
  #     dimensions; visually verified after rendering, not derived exactly.
  suppressMessages(library(png))
  PANEL_W_IN <- 4.2; PANEL_H_IN <- 6.0  # rough per-panel physical size estimate
  # corner = "bottomleft"/"topright" anchors to that plot corner; corner =
  # "topleft_at" anchors the image's top-left corner to a specific data
  # point instead (anchor = c(x, y)), for placement next to a specific
  # contour label rather than a plot corner.
  add_silhouette <- function(p, img_path, corner, xrange, yrange, height_frac = 0.35, anchor = NULL) {
    img <- png::readPNG(img_path)
    img_h_px <- dim(img)[1]; img_w_px <- dim(img)[2]
    h_data <- height_frac * diff(yrange)
    h_in   <- height_frac * PANEL_H_IN
    w_in   <- h_in * (img_w_px / img_h_px)
    w_data <- w_in * (diff(xrange) / PANEL_W_IN)
    if (corner == "bottomleft") {
      xmin <- xrange[1]; xmax <- xmin + w_data
      ymin <- yrange[1]; ymax <- ymin + h_data
    } else if (corner == "topright") {
      xmax <- xrange[2]; xmin <- xmax - w_data
      ymax <- yrange[2]; ymin <- ymax - h_data
    } else {  # topleft_at
      xmin <- anchor[1]; xmax <- xmin + w_data
      ymax <- anchor[2]; ymin <- ymax - h_data
    }
    p + annotation_custom(grid::rasterGrob(img, interpolate = TRUE),
                           xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  }

  panel_a_xrange <- range(hw_panel_a_data$fog_modern, na.rm = TRUE)
  panel_a_yrange <- range(hw_panel_a_data$BIO18_modern, na.rm = TRUE)
  panel_a <- panel_a %>%
    add_silhouette("/private/groups/meyerlab/esteiger/pheno/tall.png", "bottomleft", panel_a_xrange, panel_a_yrange, height_frac = 0.46) %>%
    add_silhouette("/private/groups/meyerlab/esteiger/pheno/wide.png", "topright", panel_a_xrange, panel_a_yrange, height_frac = 0.46)

  # NOTE: global z_vars has since been overwritten by the historic section
  # below (it ends as c("Fog_z","BIO18_z","BIO15_z")) — modern's own set is
  # hardcoded here rather than reused, to match panel_b's actual model.
  ray_prop_panel_data <- modern %>% drop_na(ray_prop, Fog_z, BIO18_z)
  panel_b_xrange <- range(ray_prop_panel_data$fog_modern, na.rm = TRUE)
  panel_b_yrange <- range(ray_prop_panel_data$BIO18_modern, na.rm = TRUE)
  panel_b <- panel_b %>%
    add_silhouette("/private/groups/meyerlab/esteiger/pheno/lower_ratio.png", "bottomleft", panel_b_xrange, panel_b_yrange, height_frac = 0.17) %>%
    add_silhouette("/private/groups/meyerlab/esteiger/pheno/higher_ratio.png", "topleft_at", panel_b_xrange, panel_b_yrange, height_frac = 0.17,
                    anchor = c(13, panel_b_yrange[2] - 0.1 * diff(panel_b_yrange)))

  combined_4panel <- (panel_a | panel_b | panel_c | panel_d) +
    plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = ".") &
    theme(plot.tag.position = c(0.02, 0.99),
          plot.tag = element_text(size = RD_TIME_FONT_PT, hjust = 0, vjust = 1),
          plot.margin = margin(t = 5, r = 8, b = 2, l = 2))
  ggsave(file.path(combined_out_dir, "phenotype.jpg"), combined_4panel,
         width = 18, height = 7, dpi = 150)
  cat(sprintf("\nCombined 4-panel figure written to %s\n",
              file.path(combined_out_dir, "phenotype.jpg")))
} else {
  warning("Not all 4 panels (a/b/c/d) were captured — combined figure skipped. Missing: ",
          paste(setdiff(c("a", "b", "c", "d"), names(PHENOTYPE_DATA_PANELS)), collapse = ", "))
}

# ============================================================
# 11b. Reordered main figure (2026-10-07 per user decision):
#        A: Height:width response surface (modern) — was B
#        B: ray_prop response surface (modern) — was C
#        C: BIO18/Kling-fog-over-time, pooled LOESS (historic) — was D
#        D: disk/ray floret counts over time (historic) — was A
# ============================================================

if (exists("panel_a") && exists("panel_b") && exists("panel_c") && exists("panel_d")) {
  combined_4panel_reordered <- (panel_a | panel_b | panel_c | panel_d) +
    plot_layout(widths = c(1.15, 1.15, 0.85, 0.85)) +
    plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = ".") &
    theme(plot.tag.position = c(0.02, 0.99),
          plot.tag = element_text(size = RD_TIME_FONT_PT, hjust = 0, vjust = 1),
          plot.margin = margin(t = 5, r = 8, b = 2, l = 2))
  ggsave(file.path(combined_out_dir, "phenotype_reordered.jpg"), combined_4panel_reordered,
         width = 18, height = 7, dpi = 150)
  cat(sprintf("Reordered main figure written to %s\n",
              file.path(combined_out_dir, "phenotype_reordered.jpg")))
} else {
  warning("Panels a/b/c/d not all captured — reordered combined figure skipped.")
}

# ============================================================
# 12. Six-panel figure (half page, 183 mm wide):
#       row 1: A plant-architecture photographs
#              B Height:width response surface (Fog + BIO18; contours 1-5)
#              C ray-floret-proportion response surface
#       row 2: D/E Kling fog trend, California, with study-region inset (E)
#              F/G trend in occurrence of summer rain, California, with
#                  study-region inset (G)
#              H disk/ray floret counts over time
#     B/C come from PHENOTYPE_FIG_PANELS (drawn in the modern section);
#     H re-themes the year-trend GLMM panel (PHENOTYPE_DATA_PANELS$d).
#     Map panels: viridis, reversed so declines are yellow and increases
#     purple; data over a light grey California polygon, a black box
#     marking the study region, inset of that region (top right) with
#     collection localities (black points), colour bar below the inset.
#       Fog trend: per-pixel least-squares slope of the annual Kling fog
#         fraction, 1950-2014, x100 (% of hours per year).
#       Summer rain: per-pixel linear-probability slope of the occurrence
#         of measurable (> 1 mm) June-August rain, 1927-2025, x100 (%
#         per decade), from env/ca_summer_precip_trend_map.R (cached tif);
#         shown unmasked by significance.
# ============================================================

suppressMessages({ library(terra); library(sf); library(jpeg) })

FIG_H_IN <- 5.8          # figure width follows from the bottom row (see layout below)
ENV_DIR  <- "/private/groups/meyerlab/esteiger/env"
# California outline. study_area/california_boundary.gpkg (from the `maps`
# package) is too coarse at the coast -- over half the collection localities
# fall outside it -- so the outline used here is the PRISM 800-m land mask
# (true coastline) within that polygon buffered by 3 km (its straight inland
# borders are fine). Cached.
CA_LAND_GPKG <- file.path(ENV_DIR, "study_area/california_land_prism.gpkg")
if (!file.exists(CA_LAND_GPKG)) {
  ca_coarse <- st_read(file.path(ENV_DIR, "study_area/california_boundary.gpkg"), quiet = TRUE)
  ca_buf    <- ca_coarse %>% st_transform(3310) %>% st_buffer(3000) %>% st_transform(4269)
  land      <- !is.na(rast(file.path(ENV_DIR, "PRISM/202006/prism_ppt_us_30s_202006.tif")))
  land      <- mask(classify(land, cbind(0, NA)), vect(ca_buf))
  st_write(st_as_sf(as.polygons(land, dissolve = TRUE)), CA_LAND_GPKG, quiet = TRUE)
}
CA_SF <- st_read(CA_LAND_GPKG, quiet = TRUE) %>% st_transform(4326)
CA_BB <- st_bbox(CA_SF)

# --- trend rasters (computed once, then cached) ---
KLING_TREND_TIF <- file.path(ENV_DIR, "figures/kling_fog_trend_pct_hours_per_year.tif")
if (!file.exists(KLING_TREND_TIF)) {
  kling <- rast(file.path(ENV_DIR, "fog/Kling/timeseries.tif"))
  yrs   <- as.integer(names(kling))
  w     <- (yrs - mean(yrs)) / sum((yrs - mean(yrs))^2)   # OLS slope = sum(w_t * y_t)
  kling_trend <- kling[[1]] * w[1]
  for (i in 2:nlyr(kling)) kling_trend <- kling_trend + kling[[i]] * w[i]
  kling_trend <- kling_trend * 100                        # fraction -> % of hours, per year
  writeRaster(kling_trend, KLING_TREND_TIF, overwrite = TRUE)
}
# Summer-rain occurrence trend, recomputed here without masking to the coarse
# polygon (env/ca_summer_precip_trend_map.R's cached tif is masked to it):
# per pixel, OLS slope of P(June-August rain > 1 mm) on year, over every
# 1927-2025 year with all three months present, x10 (per decade) x100 (%).
RAIN_TREND_TIF <- file.path(ENV_DIR, "figures/ca_summer_rain_occurrence_trend_pct_per_decade_unmasked.tif")
if (!file.exists(RAIN_TREND_TIF)) {
  prism_f <- function(y, m) file.path(ENV_DIR, "PRISM", sprintf("%d%02d", y, m), sprintf("prism_ppt_us_30s_%d%02d.tif", y, m))
  ryrs <- Filter(function(y) all(file.exists(prism_f(y, 6:8))), 1927:2025)
  w    <- (ryrs - mean(ryrs)) / sum((ryrs - mean(ryrs))^2)
  rain_acc <- NULL
  for (i in seq_along(ryrs)) {
    jja <- Reduce(`+`, lapply(prism_f(ryrs[i], 6:8), rast))
    term <- (jja > 1) * w[i]
    rain_acc <- if (is.null(rain_acc)) term else rain_acc + term
  }
  writeRaster(rain_acc * 10 * 100, RAIN_TREND_TIF, overwrite = TRUE)
  cat(sprintf("Summer-rain occurrence trend computed from %d years (%d-%d)\n", length(ryrs), min(ryrs), max(ryrs)))
}
fog_trend  <- mask(rast(KLING_TREND_TIF), vect(CA_SF))
rain_trend <- mask(rast(RAIN_TREND_TIF), vect(st_transform(CA_SF, 4269)))
# clip to the fog raster's footprint (the Kling coastal band), so D and F
# show the same area
rain_trend <- mask(rain_trend, resample(fog_trend, rain_trend, method = "near"))

# --- collection localities and the study-region box ---
fig_pts <- read.csv("/private/groups/meyerlab/esteiger/pheno/env_wi_range.csv")[, c("Lat", "Long")] %>%
  distinct() %>% filter(!is.na(Lat), !is.na(Long))
# the project's study box (bbox_extent = ext(-120.7, -120.1, 34.4, 34.9))
STUDY_BOX <- c(xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9)
study_box_sf <- st_as_sfc(st_bbox(STUDY_BOX, crs = st_crs(4326)))

MAP_EAST_PAD <- 1.4   # degrees of empty frame east of California (inset + colour bar)
MAP_XLIM <- c(CA_BB[["xmin"]], CA_BB[["xmax"]] + MAP_EAST_PAD)
MAP_YLIM <- c(CA_BB[["ymin"]], CA_BB[["ymax"]] + 0.75)   # empty band above California for the tag
MAP_ASPECT <- diff(MAP_XLIM) * cos(mean(MAP_YLIM) * pi / 180) / diff(MAP_YLIM)   # width / height
MAP_H_IN   <- 2.56                                    # height of the map panels (bottom row)
MAP_W_IN   <- MAP_H_IN * MAP_ASPECT                   # width of D and of F
# inset slot (fractions of the map panel) and where the fixed-aspect inset's
# right edge actually falls inside it, so the colour bar can align with it
INSET_SLOT   <- c(left = 0.40, bottom = 0.60, right = 1, top = 1)
INSET_ASPECT <- (0.6 * cos(34.65 * pi / 180)) / 0.5   # study box: 0.6 deg lon x 0.5 deg lat
inset_w_in   <- min((INSET_SLOT[["right"]] - INSET_SLOT[["left"]]) * MAP_W_IN,
                    (INSET_SLOT[["top"]] - INSET_SLOT[["bottom"]]) * MAP_H_IN * INSET_ASPECT)
INSET_RIGHT  <- 1 - ((INSET_SLOT[["right"]] - INSET_SLOT[["left"]]) * MAP_W_IN - inset_w_in) / 2 / MAP_W_IN
raster_df <- function(r) { d <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(d)[3] <- "value"; d }

map_panel <- function(r, legend_title, tag_main, tag_inset) {
  lims  <- as.numeric(quantile(values(r), c(0.01, 0.99), na.rm = TRUE))
  fill  <- scale_fill_viridis_c(direction = -1, limits = lims, oob = scales::squish, name = legend_title)
  main <- ggplot() +
    geom_sf(data = CA_SF, fill = "grey92", color = "grey55", linewidth = 0.2) +
    geom_raster(data = raster_df(aggregate(r, 2, na.rm = TRUE)), aes(x = x, y = y, fill = value)) +
    geom_sf(data = study_box_sf, fill = NA, color = "black", linewidth = 0.35) +
    fill +
    guides(fill = guide_colourbar(barwidth = unit(5, "pt"), barheight = unit(48, "pt"),
                                  title.position = "top")) +
    # extend the frame east of California so the inset and colour bar sit in
    # empty space at the right of the panel rather than over the data
    coord_sf(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE, datum = NA) +
    labs(tag = tag_main) +
    theme_void(base_size = FIG_PT) +
    # colour bar directly below the inset, right-aligned with it
    theme(legend.position = "inside", legend.position.inside = c(INSET_RIGHT, 0.57),
          legend.justification = c(1, 1),                     # right edge level with the inset's
          legend.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = FIG_PT - 0.5, hjust = 0),
          legend.text  = element_text(size = FIG_PT - 1),
          plot.tag = element_text(size = FIG_PT + 2, hjust = 0, vjust = 1),
          plot.tag.position = c(0, 1),
          plot.margin  = margin(1, 1, -6, 1))       # bottom of California level with H's x-axis title
  inset <- ggplot() +
    geom_raster(data = raster_df(crop(r, ext(STUDY_BOX))), aes(x = x, y = y, fill = value)) +
    geom_point(data = fig_pts, aes(x = Long, y = Lat), shape = 21, fill = "white", color = "black",
               stroke = 0.3, size = 0.9) +                      # as in panels B/C
    fill + guides(fill = "none") +
    coord_sf(xlim = STUDY_BOX[c("xmin", "xmax")], ylim = STUDY_BOX[c("ymin", "ymax")],
             expand = FALSE, datum = NA, crs = 4326, default_crs = 4326) +
    labs(tag = tag_inset) +
    theme_void(base_size = FIG_PT) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
          panel.background = element_rect(fill = "white", color = NA),
          plot.tag = element_text(size = FIG_PT + 2, hjust = 1, vjust = 1),
          plot.tag.position = c(-0.07, 1))          # just left of the inset
  main + inset_element(inset, left = INSET_SLOT[["left"]], bottom = INSET_SLOT[["bottom"]],
                       right = INSET_SLOT[["right"]], top = INSET_SLOT[["top"]], align_to = "panel",
                       ignore_tag = TRUE)
}

panel_D <- map_panel(fog_trend,  "Fog trend\n(% hrs/yr)", "D", "E")
panel_F <- map_panel(rain_trend, "Summer rain trend\n(% pts/decade)", "F", "G")

# --- A: photographs ---
arch <- jpeg::readJPEG("/private/groups/meyerlab/esteiger/pheno/modern_data_output/plant_architecture.jpg")
arch <- arch[seq(1, dim(arch)[1], 4), seq(1, dim(arch)[2], 4), ]
panel_A <- ggplot() +
  annotation_raster(arch, xmin = 0, xmax = dim(arch)[2], ymin = 0, ymax = dim(arch)[1]) +
  scale_x_continuous(limits = c(0, dim(arch)[2]), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, dim(arch)[1]), expand = c(0, 0)) +
  coord_fixed() + labs(tag = "A") + theme_void() + theme(plot.margin = margin(1, 1, 1, 1))

# --- B, C: response surfaces with phenotype-extreme silhouettes ---
fig_silhouette <- function(p, img_path, corner, xr, yr, height_frac, anchor = NULL,
                           panel_w_in = 1.9, panel_h_in = 2.1) {
  img <- png::readPNG(img_path)
  h_data <- height_frac * diff(yr)
  w_data <- height_frac * panel_h_in * (dim(img)[2] / dim(img)[1]) * diff(xr) / panel_w_in
  if (corner == "bottomleft") { x0 <- xr[1]; y0 <- yr[1] }
  else if (corner == "topright") { x0 <- xr[2] - w_data; y0 <- yr[2] - h_data }
  else { x0 <- anchor[1]; y0 <- anchor[2] - h_data }
  p + annotation_custom(grid::rasterGrob(img, interpolate = TRUE),
                        xmin = x0, xmax = x0 + w_data, ymin = y0, ymax = y0 + h_data)
}
SIL <- "/private/groups/meyerlab/esteiger/pheno"
xrB <- range(hw_panel_a_data$fog_modern, na.rm = TRUE); yrB <- range(hw_panel_a_data$BIO18_modern, na.rm = TRUE)
panel_B <- PHENOTYPE_FIG_PANELS$B %>%
  fig_silhouette(file.path(SIL, "tall.png"), "bottomleft", xrB, yrB, 0.55) %>%
  fig_silhouette(file.path(SIL, "wide.png"), "topright",   xrB, yrB, 0.55)
cd  <- PHENOTYPE_FIG_PANELS$C_data
xrC <- range(cd$fog_modern, na.rm = TRUE); yrC <- range(cd$BIO18_modern, na.rm = TRUE)
panel_C <- PHENOTYPE_FIG_PANELS$C %>%
  fig_silhouette(file.path(SIL, "lower_ratio.png"),  "bottomleft", xrC, yrC, 0.21) %>%
  fig_silhouette(file.path(SIL, "higher_ratio.png"), "topleft_at", xrC, yrC, 0.21,
                 anchor = c(13, yrC[2] - 0.1 * diff(yrC)))
panel_B <- panel_B + labs(tag = "B"); panel_C <- panel_C + labs(tag = "C")

# --- H: floret counts over time, at figure scale ---
panel_H <- PHENOTYPE_DATA_PANELS$d
panel_H$layers[[1]]$aes_params$size <- 0.9
panel_H <- panel_H + labs(tag = "H") + FIG_THEME +
  theme(axis.title = element_text(size = FIG_PT), axis.text = element_text(size = FIG_PT - 0.5),
        legend.text = element_text(size = FIG_PT - 0.5), legend.key.size = unit(8, "pt"),
        legend.margin = margin(1, 2, 1, 2), plot.margin = margin(1, 2, 1, 3))

# tag style set per panel (not with `&`, which would also override the
# insets' own tag placement)
TAG_THEME <- theme(plot.tag = element_text(size = FIG_PT + 2, hjust = 0, vjust = 1), plot.tag.position = c(0, 1))
# Layout in absolute units so every panel fills its cell:
#   bottom row: D and F exactly as wide as their (fixed-aspect) maps,
#     side by side with no gap; H takes the rest, wrapped as a whole plot so its
#     whole plot, axis title included, spans the maps' height.
#   top row: A exactly as wide as the photograph at the row's height;
#     B and C share the remaining width.
H_W_IN    <- 2.25                                     # width of H
FIG_W_IN  <- 2 * MAP_W_IN + H_W_IN                    # bottom row sets the figure width; top row
                                                      # (A at the photo's width, B/C share the rest) matches it
TOP_H_IN  <- FIG_H_IN - MAP_H_IN - 0.12               # full top-row height (A spans B/C's axis titles too)
ARCH_W_IN <- (TOP_H_IN - 0.32) * dim(arch)[2] / dim(arch)[1]   # 0.32 in: tag line + row margins
# A wrapped as a whole plot so the photograph fills the row's full height,
# its bottom level with B/C's x-axis titles
top_row <- (wrap_elements(full = panel_A + TAG_THEME) | panel_B + TAG_THEME | panel_C + TAG_THEME) +
  plot_layout(widths = unit(c(ARCH_W_IN, 1, 1), c("in", "null", "null")))
bottom_row <- (panel_D | panel_F | wrap_elements(full = panel_H + TAG_THEME)) +
  plot_layout(widths = unit(c(MAP_W_IN, MAP_W_IN, 1), c("in", "in", "null")))
fig6 <- top_row / bottom_row + plot_layout(heights = unit(c(1, MAP_H_IN), c("null", "in")))
fig_out <- file.path(combined_out_dir, "phenotype_figure")
ggsave(paste0(fig_out, ".jpg"), fig6, width = FIG_W_IN, height = FIG_H_IN, dpi = 300)
cat(sprintf("Six-panel figure written to %s.jpg (%.1f x %.1f in)\n", fig_out, FIG_W_IN, FIG_H_IN))
save(PHENOTYPE_FIG_PANELS, PHENOTYPE_DATA_PANELS, hw_panel_a_data,
     file = file.path(combined_out_dir, "phenotype_figure_panels.RData"))
