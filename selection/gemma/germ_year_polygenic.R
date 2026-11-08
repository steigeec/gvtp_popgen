#!/usr/bin/env Rscript
# Germination-year tests for a variable's 40000_3 outlier peaks in Strauss: per-peak ANOVA and
# regression of favorable-allele dosage on germination-year values, focal peak plots, and a
# polygenic score across peaks.
# Usage: Rscript germ_year_polygenic.R <var_name> <long_term_col> <germ_yr_col>

BASE_DIR <- "/private/groups/meyerlab/esteiger"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript germ_year_polygenic.R <var_name> <long_term_col> <germ_yr_col>")
}

var_name       <- args[1]
long_term_col  <- args[2]
germ_yr_col    <- args[3]

.libPaths(file.path(BASE_DIR, "software/R_libs"))
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)

env_labels <- c(
  fog     = "Fog (days per month)",
  BIO12   = "Annual precipitation (mm)",
  BIO15   = "Precipitation seasonality",
  BIO18   = "Precipitation in the\nwarmest quarter (mm)",
  BCM_30y = "Climate water deficit (mm)"
)

year_fills <- c("2020" = "white", "2021" = "grey30", "2022" = "grey70")

base_theme <- theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        strip.text = element_text(size = 11))

MIN_SNPS_POLYGENIC <- 30
N_FOCAL_PEAKS      <- 5

base_dir     <- file.path(BASE_DIR, "selection/gemma")
env_file     <- file.path(BASE_DIR, "env_var.txt")
scaffold_map <- file.path(base_dir, "scaffold_snp_map.txt")
out_dir      <- file.path(base_dir, var_name)

report_path <- file.path(out_dir, paste0("germ_year_polygenic_", var_name, "_report.txt"))
sink(report_path, split = TRUE)

# ── Read outlier SNPs ────────────────────────────────────────────────────────
outlier_file <- file.path(out_dir, paste0("outliers_for_annot_40000_3_", var_name, ".txt"))
if (!file.exists(outlier_file)) {
  outlier_file <- file.path(out_dir, paste0("unfiltered_outliers_for_annot_", var_name, ".txt"))
}
if (!file.exists(outlier_file)) {
  cat("No outlier file found for", var_name, "-- nothing to do.\n")
  sink()
  quit(save = "no", status = 0)
}

outliers <- fread(outlier_file, header = FALSE)
snp_info <- data.frame(
  snp = as.character(outliers$V2),
  chr = outliers$V1,
  pos = outliers$V3,
  stringsAsFactors = FALSE
)
snp_info <- snp_info[!duplicated(snp_info$snp), ]
snp_info <- snp_info[order(snp_info$chr, snp_info$pos), ]
cat("Total outlier SNPs:", nrow(snp_info), "\n")

# ── Identify peaks (merge SNPs within 40kb on same chromosome) ───────────────
snp_info$peak <- 1L
if (nrow(snp_info) > 1) {
  for (i in 2:nrow(snp_info)) {
    if (snp_info$chr[i] != snp_info$chr[i - 1] ||
        snp_info$pos[i] - snp_info$pos[i - 1] > 40000) {
      snp_info$peak[i] <- snp_info$peak[i - 1] + 1L
    } else {
      snp_info$peak[i] <- snp_info$peak[i - 1]
    }
  }
}

peak_sizes <- table(snp_info$peak)
peaks_3plus <- as.integer(names(peak_sizes[peak_sizes >= 3]))
snp_info_peaks <- snp_info[snp_info$peak %in% peaks_3plus, ]
cat("Peaks with >= 3 SNPs:", length(peaks_3plus), "\n")

# ── Map numeric chromosome IDs to DEINVI_* names ────────────────────────────
smap <- fread(scaffold_map, select = c(1, 3), col.names = c("chr_num", "chr_name"))
smap <- unique(smap)
chr_lookup <- smap$chr_name
names(chr_lookup) <- as.character(smap$chr_num)

# ── Read dosages ─────────────────────────────────────────────────────────────
dosage_file <- file.path(out_dir, paste0("dosage_", var_name, ".txt"))
if (!file.exists(dosage_file)) {
  stop("Dosage file not found: ", dosage_file,
       "\nRun the DS extraction step in Launch_find_outlier_annotate.sh first.")
}

cat("Reading DS dosages from", dosage_file, "\n")
ds <- fread(dosage_file, header = TRUE)
ds_snps <- as.character(ds$snp)
ds_samples <- colnames(ds)[-1]
dosage <- as.matrix(ds[, -1])
rownames(dosage) <- ds_snps
dosage <- t(dosage)

cat("Dosage matrix:", nrow(dosage), "individuals x", ncol(dosage), "SNPs\n")

# ── Read environmental data ──────────────────────────────────────────────────
env <- fread(env_file, header = TRUE, sep = "\t")
env$SampleID <- as.character(env$SampleID)
# 13 individuals have untrustworthy coordinates (Location="Strauss" but
# Clusters assigned to Point_Conception/VAFB_Lion_Head); excluded throughout,
# including from allele polarization
BAD_COORD_IDS <- c("4548", "4549", "4550", "4551", "4552", "4553", "4554",
                   "4555", "4556", "4557", "4559", "4561", "4562")
env <- env[!env$SampleID %in% BAD_COORD_IDS, ]

env_match <- match(rownames(dosage), env$SampleID)
matched <- !is.na(env_match)
dosage <- dosage[matched, , drop = FALSE]
env_matched <- env[env_match[matched], ]

cat("Samples matched to env data:", nrow(dosage), "\n")

long_term_env <- env_matched[[long_term_col]]
germ_env      <- env_matched[[germ_yr_col]]
sample_loc    <- as.character(env_matched$Location)
sample_year   <- env_matched$Year
sample_ids    <- as.character(env_matched$SampleID)

# ── Polarize alleles (vectorized regression controlling for coverage) ────────
coverage <- as.numeric(env_matched$Cov)

cov_mat <- cbind(1, coverage)
hat <- cov_mat %*% solve(crossprod(cov_mat), t(cov_mat))
env_resid <- as.numeric(long_term_env - hat %*% long_term_env)
dos_resid <- dosage - hat %*% dosage
betas <- as.numeric(crossprod(dos_resid, env_resid) / colSums(dos_resid^2))

flip <- which(betas < 0)
if (length(flip) > 0) {
  dosage[, flip] <- 2 - dosage[, flip]
}
cat("SNPs flipped to polarize:", length(flip), "of", ncol(dosage), "\n")

# ── Subset to Strauss + add germination year ─────────────────────────────────
strauss_idx <- which(sample_loc == "Strauss")

df <- data.frame(
  sample_id = sample_ids[strauss_idx],
  germ_env  = germ_env[strauss_idx],
  Year      = factor(sample_year[strauss_idx]),
  stringsAsFactors = FALSE
)

cat("Strauss samples:", nrow(df), "\n")
cat("Samples per germination year:\n")
print(table(df$Year))

# ── Per-peak analysis (ANOVA + regression) ───────────────────────────────────
strauss_dosage <- dosage[strauss_idx, , drop = FALSE]

year_env <- df %>%
  group_by(Year) %>%
  summarise(env_val = mean(germ_env), .groups = "drop")

peak_results <- list()
for (pk in peaks_3plus) {
  pk_snps <- snp_info_peaks$snp[snp_info_peaks$peak == pk]
  pk_chr  <- snp_info_peaks$chr[snp_info_peaks$peak == pk][1]
  pk_start <- min(snp_info_peaks$pos[snp_info_peaks$peak == pk])
  pk_end   <- max(snp_info_peaks$pos[snp_info_peaks$peak == pk])

  pk_cols <- which(colnames(strauss_dosage) %in% pk_snps)
  if (length(pk_cols) < 3) next

  pk_score <- rowMeans(strauss_dosage[, pk_cols, drop = FALSE], na.rm = TRUE)

  pk_aov <- aov(pk_score ~ df$Year)
  pk_p_anova <- summary(pk_aov)[[1]]["df$Year", "Pr(>F)"]

  pk_lm <- lm(pk_score ~ df$germ_env)
  pk_p_reg <- summary(pk_lm)$coefficients["df$germ_env", "Pr(>|t|)"]
  pk_beta  <- coef(pk_lm)["df$germ_env"]

  peak_results[[as.character(pk)]] <- data.frame(
    peak    = pk,
    chr     = pk_chr,
    start   = pk_start,
    end     = pk_end,
    n_snps  = length(pk_cols),
    p_anova = pk_p_anova,
    p_reg   = pk_p_reg,
    beta    = pk_beta
  )
}

if (length(peak_results) == 0) {
  cat("No valid peaks found -- skipping.\n")
  sink()
  quit(save = "no", status = 0)
}

peak_summary <- do.call(rbind, peak_results)

cat("\n-- Peak summary (top 10 by ANOVA p-value) --\n")
print(head(peak_summary[order(peak_summary$p_anova), ], 10))

cat("\n-- Peak summary (top 10 by regression p-value) --\n")
print(head(peak_summary[order(peak_summary$p_reg), ], 10))

# ── Helper: build focal peak plot ────────────────────────────────────────────
make_focal_plot <- function(focal_pks, p_col, label_prefix, filename) {
  env_ylab <- if (var_name %in% names(env_labels)) env_labels[var_name] else var_name

  p_env <- ggplot(df, aes(x = Year, y = germ_env, fill = Year)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
    scale_fill_manual(values = year_fills) +
    labs(x = "Germination Year", y = env_ylab) +
    base_theme

  focal_data <- list()
  for (i in 1:nrow(focal_pks)) {
    pk <- focal_pks$peak[i]
    pk_snps <- snp_info_peaks$snp[snp_info_peaks$peak == pk]
    pk_cols <- which(colnames(strauss_dosage) %in% pk_snps)
    pk_score <- rowMeans(strauss_dosage[, pk_cols, drop = FALSE], na.rm = TRUE)

    chr_name <- chr_lookup[as.character(focal_pks$chr[i])]
    label <- paste0(chr_name,
                    " (", focal_pks$n_snps[i], " SNPs)\n",
                    "p = ", signif(focal_pks[[p_col]][i], 3))

    focal_data[[i]] <- data.frame(
      score      = pk_score,
      Year       = df$Year,
      peak_label = label,
      stringsAsFactors = FALSE
    )
  }

  focal_df <- do.call(rbind, focal_data)
  focal_df$peak_label <- factor(focal_df$peak_label,
                                levels = sapply(focal_data, function(x) x$peak_label[1]))

  n_pk <- nrow(focal_pks)
  n_c <- min(n_pk + 1, 4)
  n_r <- ceiling((n_pk + 1) / n_c)

  p_peaks <- ggplot(focal_df, aes(x = Year, y = score, fill = Year)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.15, height = 0, alpha = 0.3, size = 0.8) +
    scale_fill_manual(values = year_fills) +
    facet_wrap(~ peak_label, scales = "free_y") +
    labs(x = "Germination Year", y = "Mean favorable allele dosage") +
    base_theme

  combined <- plot_grid(p_env, p_peaks,
                        rel_widths = c(1, min(n_pk, n_c - 1)))

  ggsave(file.path(out_dir, filename),
         combined,
         width = 4 * n_c, height = 4 * n_r)

  cat(label_prefix, "focal plot saved:", filename, "\n")
}

# ── ANOVA-ranked focal peaks ────────────────────────────────────────────────
focal_anova <- head(peak_summary[order(peak_summary$p_anova), ], N_FOCAL_PEAKS)
cat("\nTop", nrow(focal_anova), "peaks by ANOVA:\n")
print(focal_anova)

cat("\n-- Tukey HSD for ANOVA-ranked peaks --\n")
for (i in 1:nrow(focal_anova)) {
  pk <- focal_anova$peak[i]
  pk_snps <- snp_info_peaks$snp[snp_info_peaks$peak == pk]
  pk_cols <- which(colnames(strauss_dosage) %in% pk_snps)
  pk_score <- rowMeans(strauss_dosage[, pk_cols, drop = FALSE], na.rm = TRUE)
  chr_name <- chr_lookup[as.character(focal_anova$chr[i])]
  cat("\n", chr_name, " (", focal_anova$n_snps[i], " SNPs)\n", sep = "")
  aov_fit <- aov(pk_score ~ df$Year)
  cat("  ANOVA p =", signif(summary(aov_fit)[[1]]["df$Year", "Pr(>F)"], 3), "\n")
  print(TukeyHSD(aov_fit))
}

make_focal_plot(focal_anova, "p_anova", "ANOVA",
                paste0("germ_year_focal_peaks_anova_", var_name, ".jpg"))

# ── Regression-ranked focal peaks ───────────────────────────────────────────
focal_reg <- head(peak_summary[order(peak_summary$p_reg), ], N_FOCAL_PEAKS)
cat("\nTop", nrow(focal_reg), "peaks by regression:\n")
print(focal_reg)

make_focal_plot(focal_reg, "p_reg", "Regression",
                paste0("germ_year_focal_peaks_regression_", var_name, ".jpg"))

# ── Polygenic score (only if enough peaks) ───────────────────────────────────
peak_rep_snps <- snp_info_peaks %>%
  group_by(peak) %>%
  slice(1) %>%
  pull(snp)

if (length(peak_rep_snps) >= MIN_SNPS_POLYGENIC) {
  cat("\n-- Polygenic score analysis (", length(peak_rep_snps),
      "peak-representative SNPs >= threshold", MIN_SNPS_POLYGENIC, ") --\n")

  rep_cols <- which(colnames(strauss_dosage) %in% peak_rep_snps)
  strauss_rep <- strauss_dosage[, rep_cols, drop = FALSE]
  n_scored   <- rowSums(!is.na(strauss_rep))
  poly_score <- rowSums(strauss_rep, na.rm = TRUE) / n_scored
  df$poly_score <- poly_score

  model <- lm(poly_score ~ germ_env, data = df)
  cat("\n-- Linear model: poly_score ~ germ_env --\n")
  print(summary(model))

  anova_res <- aov(poly_score ~ Year, data = df)
  cat("\n-- ANOVA: poly_score ~ Year --\n")
  print(summary(anova_res))

  tukey <- TukeyHSD(anova_res)
  cat("\n-- Tukey HSD --\n")
  print(tukey)

  # ── Directed test: is poly_score higher in the highest-germ_env year? ───────
  # One-sided; the high year is the year with the highest mean germ_env
  year_means <- df %>% group_by(Year) %>% summarise(mean_germ_env = mean(germ_env), .groups = "drop")
  cat("\n-- Mean", germ_yr_col, "by Year --\n")
  print(year_means)
  high_year <- as.character(year_means$Year[which.max(year_means$mean_germ_env)])
  cat("Highest", germ_yr_col, "year:", high_year, "\n")

  cat("\n-- One-sided t-test: poly_score in", high_year, "> poly_score in pooled other years --\n")
  print(t.test(df$poly_score[df$Year == high_year], df$poly_score[df$Year != high_year],
               alternative = "greater"))

  # Violin plot by germination year
  p1 <- ggplot(df, aes(x = Year, y = poly_score, fill = Year)) +
    geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75), color = "black") +
    geom_jitter(width = 0.15, alpha = 0.2, size = 0.5) +
    scale_fill_manual(values = year_fills) +
    labs(x = "Germination Year",
         y = "Mean favorable allele dosage") +
    base_theme

  ggsave(file.path(out_dir, paste0("germ_year_polygenic_", var_name, ".jpg")),
         p1, width = 6, height = 5)

  # Scatterplot colored by year
  r2   <- summary(model)$r.squared
  pval <- summary(model)$coefficients[2, 4]

  p2 <- ggplot(df, aes(x = germ_env, y = poly_score, fill = Year)) +
    geom_jitter(shape = 21, color = "black", alpha = 1,
                width = diff(range(df$germ_env)) * 0.02, size = 1.5, stroke = 0.3) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
    scale_fill_manual(values = year_fills) +
    labs(x = paste0(var_name, " (germination year)"),
         y = "Mean favorable allele dosage") +
    annotate("text",
             x = min(df$germ_env), y = max(df$poly_score),
             label = paste0("R2 = ", round(r2, 4), "\np = ",
                            format.pval(pval, digits = 3)),
             hjust = 0, vjust = 1, size = 5) +
    base_theme

  ggsave(file.path(out_dir, paste0("germ_year_scatter_", var_name, ".jpg")),
         p2, width = 6, height = 5)

  cat("Polygenic score plots saved.\n")
} else {
  cat("\nSkipping polygenic score (", length(peak_rep_snps),
      "peak-representative SNPs < threshold", MIN_SNPS_POLYGENIC, ").\n")
}

cat("\nDone.\n")
sink()
