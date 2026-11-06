# fire_turbine_analyses.R
#
# Standalone home for this project's three iNat-dependent "secondary"
# analyses, split out of phenotype_data.R (2026-10-06, per user
# request) so that script's main pipeline/combined figure no longer
# carries fire- or turbine-specific content, and so the one-time
# ray-count validation test isn't duplicated across modern/historic
# sections there:
#   1. Ray-floret count validation: paired Wilcoxon test justifying
#      treating dissection counts and iNat-photo counts as
#      interchangeable (phenotype_data.R's own ray-count fusion still
#      happens there, independently, on both the modern and historic
#      datasets — this is just the one-time justification for it).
#   2. Fire history effect on floret counts/proportion (historic data;
#      years since the most recent fire before collection, from
#      extract_fire_data.R / fire_data.txt).
#   3. Wind turbine effect on ray floret count at Strauss (iNat-
#      preferred ray-count fusion, difference-in-differences GLMM vs.
#      two control localities).
#
# Fire and turbine each rebuild their own copy of the phenotype dataset
# they need (historic-style for fire, a 3-locality iNat-preferred fusion
# for turbine) rather than reusing phenotype_data.R's `modern`/
# `historic` objects — duplicated here rather than sourced, since
# phenotype_data.R is a long script that does a lot of unrelated
# modeling/plotting neither analysis needs.
#
# iNat data (inatrays.csv) still contributes to phenotype_data.R
# independent of this split — ray-count fusion (the actual merge, as
# opposed to its validation) and locality backfill still happen there,
# in both the modern and historic sections.
#
# Outputs:
#   historic_data_output/fire_effect_report.csv    — GLMM (burned-only, primary)
#                                                      + sentinel-coded GLMM and
#                                                      Wilcoxon (supplementary)
#   historic_data_output/fire_effect_boxplots.png  — burned vs. unburned, by phenotype
#   historic_data_output/fire_effect_scatter.png   — years since fire vs. phenotype
#   modern_data_output/turbine_report.csv          — DiD GLMM coefficients
#   modern_data_output/turbine_ray_before_after_2020.png
#   modern_data_output/turbine_map.png

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(glmmTMB)
  library(ggplot2)
})

# ============================================================
# Paths
# ============================================================

env_file      <- "/private/groups/meyerlab/esteiger/pheno/emma_gee_results_final.csv"
pheno_file    <- "/private/groups/meyerlab/esteiger/pheno/gtinatandvars.csv"
inatrays_file <- "/private/groups/meyerlab/esteiger/pheno/inatrays.csv"
fire_file     <- "/private/groups/meyerlab/esteiger/pheno/fire_data.txt"
out_dir       <- "/private/groups/meyerlab/esteiger/pheno/historic_data_output"
turbine_out_dir <- "/private/groups/meyerlab/esteiger/pheno/modern_data_output"
dir.create(out_dir, showWarnings = FALSE)
dir.create(turbine_out_dir, showWarnings = FALSE)

SIG_ALPHA          <- 0.05
PLOT_BASE_SIZE     <- 20
CONTOUR_LABEL_SIZE <- 7
CONTOUR_LABEL_PT   <- CONTOUR_LABEL_SIZE * ggplot2::.pt

normalize_locality <- function(x) {
  x <- gsub("_", "", x)
  as.character(ifelse(is.na(x) | x == "", NA_character_, x))
}

# ============================================================
# 1. Ray-floret count validation: paired Wilcoxon test
#
#    The undergrad's stattests.R set up (but didn't capture output for)
#    a paired comparison of dissection-based Count_ray_flourets vs.
#    iNat_ray_count for specimens measured both ways, to justify using
#    one as a stand-in for the other. Run here, once, on the full
#    phenotype file (not year-restricted — matches the scope their own
#    script used) as a paired Wilcoxon signed-rank test (Shapiro-Wilk on
#    the differences rejects normality, matching their own "#not
#    normal" comment, so a paired t-test isn't appropriate): if not
#    significant, phenotype_data.R treats dissection and iNat ray counts
#    as interchangeable in its own, separate ray-count fusion (dissection
#    where available, iNat count otherwise — or the reverse precedence
#    used below for the turbine analysis).
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
            ") — fusing them in phenotype_data.R is NOT well-justified by this test.")
  } else {
    cat("  Not significant — dissection and iNat ray counts are treated as interchangeable.\n\n")
  }
}

# ============================================================
# 2. Rebuild the historic phenotype dataset (identical logic to
#    phenotype_data.R's historic Sections 1-2 — see there for the full
#    rationale; not reproduced here)
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

historic <- inner_join(historic, env, by = "SampleID")

inat <- read.csv(inatrays_file, stringsAsFactors = FALSE)
inat$SampleID <- trimws(as.character(inat$inat_Specimen))
inat_expansion <- inat %>%
  filter(!is.na(iNat_ray_count), !is.na(latitude), !is.na(longitude),
         !(SampleID %in% historic$SampleID)) %>%
  transmute(SampleID, year = Year, Locality = normalize_locality(Locality),
            Count_disk_flourets = NA_real_, Count_ray_flourets = iNat_ray_count) %>%
  distinct() %>%
  inner_join(env, by = "SampleID")

historic <- bind_rows(historic, inat_expansion)

na_loc <- is.na(historic$Locality)
if (any(na_loc)) {
  inat_loc_by_id <- setNames(normalize_locality(inat$Locality), inat$SampleID)
  historic$Locality[na_loc] <- unname(inat_loc_by_id[historic$SampleID[na_loc]])
}

historic$Locality <- droplevels(as.factor(historic$Locality))
historic$ray_prop <- with(historic, Count_ray_flourets / (Count_disk_flourets + Count_ray_flourets))
historic$Count_total_florets <- with(historic, Count_disk_flourets + Count_ray_flourets)

cat(sprintf("Historic phenotype dataset rebuilt: %d records\n", nrow(historic)))

# ============================================================
# 3. GLMM helpers (identical to phenotype_data.R's Section 4 — see
#    there for the full rationale; only fit_one/get_coefs are needed
#    here, since this script fits fixed formulas rather than ranking a
#    candidate set)
# ============================================================

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

# ============================================================
# 4. Fire history effects on floret counts/proportion (historic)
#
#    Adapted from the undergrad's fire.R (extraction-only portion
#    re-implemented in extract_fire_data.R, since fire.R mixed
#    extraction/GLMs/map-making in one file with hardcoded Google
#    Drive paths and relied on `merged`/`camap` objects never defined
#    anywhere in that file). fire_data.txt has, per specimen: years
#    since the most recent fire before collection (years_since_fire,
#    NA if never burned — not yet sentinel-coded), burned Y/N, that
#    fire's acreage/cause, and a count of fires before collection.
#
#    fire.R modeled "never burned" as years_since_fire = -250 (a
#    deliberately extreme "long enough ago not to matter" sentinel,
#    not a real measurement) so glmer() wouldn't just drop those rows
#    via NA. Kept here for continuity with the original approach —
#    years_since_fire is signed so more-recent fires are smaller
#    positive numbers (2026-09-14 convention: collection year minus
#    fire year), so the "never burned" sentinel is +250 to preserve the
#    same "long enough ago not to matter" meaning. It means the GLM's
#    "years since fire" coefficient is really being estimated off the
#    burned/unburned CONTRAST for most of its power (few specimens have
#    a real, non-sentinel value) — the Wilcoxon burned-vs-unburned
#    tests below are the more literal read of the same question.
# ============================================================

cat("\n=== Fire history effects on floret counts/proportion ===\n\n")

fire_data <- read.csv(fire_file, stringsAsFactors = FALSE)
fire_data$SampleID <- trimws(as.character(fire_data$SampleID))
fire_data <- fire_data %>%
  transmute(SampleID, years_since_fire,
            years_since_fire_sentinel = ifelse(is.na(years_since_fire), 250, years_since_fire),
            burned = factor(burned, levels = c("NO", "YES")))

historic_fire <- historic %>% inner_join(fire_data, by = "SampleID")

# Only localities with any burned specimens carry information about a
# fire effect at all — confirmed here from the data itself (matches
# fire.R's own restriction to VAFBArguello/CDFW, but re-derived rather
# than assumed): every other locality has zero burned specimens, so
# they contribute nothing but noise to a burned-vs-unburned contrast.
burned_locs <- historic_fire %>% filter(burned == "YES") %>% pull(Locality) %>% unique() %>% as.character()
cat("Localities with any burned-before-collection specimens:", paste(burned_locs, collapse = ", "), "\n")
fire_sub <- historic_fire %>% filter(as.character(Locality) %in% burned_locs)
fire_sub$Locality <- droplevels(as.factor(fire_sub$Locality))
cat(sprintf("Fire-analysis specimens (those localities only): %d (%d burned, %d unburned)\n\n",
            nrow(fire_sub), sum(fire_sub$burned == "YES"), sum(fire_sub$burned == "NO")))

fire_report_rows <- list()

# Primary fire test (2026-09-14 per user decision): years_since_fire on
# real values, restricted to burned specimens only — the sentinel-coded
# version (below, kept as a supplementary check) compresses this effect
# toward zero, since one model has to fit a slope spanning from +250
# (the "never burned" placeholder) down to real values around 8-24.
# Checked: at this n (14-20, split across only 3 localities, e.g.
# 10/3/1 for ray_prop) the (1|Locality) term is typically singular
# (RE variance ~0) — the "random effect" isn't doing real adjustment
# work here, just a fixed-effects regression on a small, unevenly-split
# sample. Reported as the primary result anyway (per-locality pooling
# is still the right model to specify even when the data can't
# distinguish it from a simple regression), but this is also why the
# scatter plot (below) drops the confidence ribbon — a CI derived from
# a singular fit would imply more precision than 14 points support.
run_fire_response_burned_only <- function(response, label, family, data) {
  d <- data %>% filter(burned == "YES") %>% drop_na(!!sym(response), years_since_fire)
  d$Locality <- droplevels(d$Locality)
  cat(sprintf("--- %s (burned-only, n=%d, %d localities) ---\n", label, nrow(d), nlevels(d$Locality)))
  m <- fit_one(paste(response, "~ years_since_fire + (1 | Locality)"), d, family)
  if (!is.null(m)) {
    cc <- get_coefs(m) %>% filter(term == "years_since_fire")
    cat(sprintf("  GLMM: years_since_fire beta=%.4f se=%.4f z=%.3f p=%.3g%s\n",
                cc$estimate, cc$std_error, cc$statistic, cc$p_value,
                if (cc$p_value < SIG_ALPHA) "  *" else ""))
    fire_report_rows[[length(fire_report_rows) + 1]] <<- data.frame(
      Phenotype = label, N = nrow(d), Test = "GLMM (years since fire, burned-only) — primary",
      Estimate_SE = sprintf("%.4f ± %.4f", cc$estimate, cc$std_error),
      statistic = round(cc$statistic, 3), p_value = signif(cc$p_value, 3))
  } else {
    cat("  GLMM failed to converge (likely too few burned specimens/localities for a random effect).\n")
  }
  cat("\n")
}

run_fire_response_burned_only("Count_disk_flourets", "Count disk florets", "poisson", fire_sub)
run_fire_response_burned_only("Count_ray_flourets", "Count ray florets", "poisson", fire_sub)
run_fire_response_burned_only("ray_prop", "Proportion ray florets", "beta", fire_sub)

cat("=== Supplementary: sentinel-coded (all specimens, never-burned = +250) + Wilcoxon burned-vs-unburned ===\n\n")

run_fire_response <- function(response, label, family, data) {
  d <- data %>% drop_na(!!sym(response))
  d$Locality <- droplevels(d$Locality)
  m <- fit_one(paste(response, "~ years_since_fire_sentinel + (1 | Locality)"), d, family)
  cat(sprintf("--- %s ---\n", label))
  if (!is.null(m)) {
    cc <- get_coefs(m) %>% filter(term == "years_since_fire_sentinel")
    cat(sprintf("  GLMM (n=%d): years_since_fire beta=%.4f se=%.4f z=%.3f p=%.3g%s\n",
                nobs(m), cc$estimate, cc$std_error, cc$statistic, cc$p_value,
                if (cc$p_value < SIG_ALPHA) "  *" else ""))
    fire_report_rows[[length(fire_report_rows) + 1]] <<- data.frame(
      Phenotype = label, N = nobs(m), Test = "GLMM (years since fire, sentinel-coded) — supplementary",
      Estimate_SE = sprintf("%.4f ± %.4f", cc$estimate, cc$std_error),
      statistic = round(cc$statistic, 3), p_value = signif(cc$p_value, 3))
  } else {
    cat("  GLMM failed to converge.\n")
  }
  wt <- tryCatch(wilcox.test(as.formula(paste(response, "~ burned")), data = d), error = function(e) NULL)
  if (!is.null(wt)) {
    cat(sprintf("  Wilcoxon (burned vs. unburned, n=%d): W=%.1f p=%.3g%s\n",
                nrow(d), wt$statistic, wt$p.value, if (wt$p.value < SIG_ALPHA) "  *" else ""))
    fire_report_rows[[length(fire_report_rows) + 1]] <<- data.frame(
      Phenotype = label, N = nrow(d), Test = "Wilcoxon (burned vs unburned) — supplementary",
      Estimate_SE = NA, statistic = round(unname(wt$statistic), 3), p_value = signif(wt$p.value, 3))
  }
  cat("\n")
}

run_fire_response("Count_disk_flourets", "Count disk florets", "poisson", fire_sub)
run_fire_response("Count_ray_flourets", "Count ray florets", "poisson", fire_sub)
run_fire_response("ray_prop", "Proportion ray florets", "beta", fire_sub)

fire_report <- bind_rows(fire_report_rows)
fire_report_path <- file.path(out_dir, "fire_effect_report.csv")
write.csv(fire_report, fire_report_path, row.names = FALSE)
cat(sprintf("Fire-effect report written to %s\n", fire_report_path))

# ============================================================
# 4. Standalone figures
# ============================================================

# Boxplot, styled to match this project's established convention (no
# title/gridlines/legend redundancy, greyscale) rather than fire.R's
# own colored (darkgoldenrod1/darkorange3) boxplots.
fire_box_data <- fire_sub %>%
  transmute(burned, `Disk florets` = Count_disk_flourets, `Ray florets` = Count_ray_flourets,
            `Ray proportion` = ray_prop) %>%
  pivot_longer(-burned, names_to = "metric", values_to = "value") %>%
  drop_na(value) %>%
  mutate(metric = factor(metric, levels = c("Disk florets", "Ray florets", "Ray proportion")))

p_fire <- ggplot(fire_box_data, aes(x = burned, y = value)) +
  geom_boxplot(fill = "grey85", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 1.5) +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = "Burned before collection", y = NULL) +
  theme_bw(base_size = PLOT_BASE_SIZE) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = CONTOUR_LABEL_PT * 0.8))
ggsave(file.path(out_dir, "fire_effect_boxplots.png"), p_fire, width = 12, height = 5, dpi = 150)
cat(sprintf("Fire-effect boxplots written to %s\n", file.path(out_dir, "fire_effect_boxplots.png")))

# Scatter + lm trend line, years since fire vs. floret morphology —
# fire.R's `one`/`two` plots (scatter + geom_smooth(method="lm")).
# Restricted to actually-burned specimens (real years_since_fire, not
# the +250 sentinel) — fire.R plotted a sentinel-inclusive version too,
# but with only ~20 burned specimens against ~250+ sentinel-coded ones,
# almost every point would stack at two extremes rather than showing a
# real relationship; the burned-only version (fire.R's "two", zoomed to
# where the real data live) is the one actually worth showing. No CI
# ribbon — at this n (14-20 specimens, 3 localities, badly unbalanced)
# the matching GLMM's (1|Locality) term is singular (RE variance ~0),
# so an interval here would imply more precision than the data support.
fire_scatter_data <- fire_sub %>%
  filter(burned == "YES") %>%
  transmute(years_since_fire, `Disk florets` = Count_disk_flourets, `Ray florets` = Count_ray_flourets,
            `Ray proportion` = ray_prop) %>%
  pivot_longer(-years_since_fire, names_to = "metric", values_to = "value") %>%
  drop_na(value) %>%
  mutate(metric = factor(metric, levels = c("Disk florets", "Ray florets", "Ray proportion")))

p_fire_scatter <- ggplot(fire_scatter_data, aes(x = years_since_fire, y = value)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.6, se = FALSE) +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = "Years since fire (most recent, before collection)", y = NULL) +
  theme_bw(base_size = PLOT_BASE_SIZE) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = CONTOUR_LABEL_PT * 0.8))
ggsave(file.path(out_dir, "fire_effect_scatter.png"), p_fire_scatter, width = 12, height = 5, dpi = 150)
cat(sprintf("Fire-effect scatter (burned specimens only) written to %s\n",
            file.path(out_dir, "fire_effect_scatter.png")))

# ============================================================
# 5. Wind turbine effect on ray floret count at Strauss
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
#        (dissection first in phenotype_data.R) — because a before/after
#        design is specifically vulnerable to the measurement method
#        itself shifting between periods; preferring one method uniformly
#        where possible keeps that consistent across the comparison. (The
#        paired-validation test in Section 1 above found no significant
#        difference between the two methods either way, so this doesn't
#        contradict that result — it's a design choice, not a correction.)
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

# iNat-preferred ray count fusion (opposite precedence from
# phenotype_data.R's own ray-count fusion — see rationale above), built
# directly from pheno + inat.
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
                     file.path(turbine_out_dir, "turbine_ray_before_after_2020.png"))

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
ggsave(file.path(turbine_out_dir, "turbine_map.png"), p_map, width = 9, height = 4.5, dpi = 150)

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
turbine_report_path <- file.path(turbine_out_dir, "turbine_report.csv")
write.csv(turbine_report, turbine_report_path, row.names = FALSE)
cat(sprintf("\nTurbine report written to %s\n", turbine_report_path))
