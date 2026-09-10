#!/usr/bin/env Rscript

# Two-panel greyscale figure for the DEINVI_00006 putative inversion peak,
# for any variable whose outliers overlap it:
#   Left  – variable's germination-year values by germination cohort
#   Right – mean favorable-allele dosage at the DEINVI_00006 peak
# Plus a 3-panel map figure: Strauss-population haplotype frequency
# (standard vs. inverted, as a black/white pie) over a Tranquillon
# Mountain contour/hillshade basemap, one panel per sampling year.
#
# Currently run for fog only — of BIO15/BIO18/fog, fog is the only variable
# whose 40000_3 outliers still include a peak inside the inversion (43-50Mb
# on DEINVI_00006); BIO18 never overlapped it, and BIO15's peak there
# (previously 66 SNPs @ 48860682-49175215) dropped out of its outlier set
# entirely as of the 2026-07-20 outlier re-call — its only chr6 outlier
# peak now sits at ~3.44Mb, unrelated to this inversion. Re-add BIO15 to
# VAR_CONFIG below if a future outlier re-call brings back a peak near
# TARGET_START.
#
# The map figure is Strauss-only: it's the one site sampled across all
# three years (2020/2021/2022); the other six sites in env_var.txt were
# only sampled in 2020, so they can't contribute to a per-year comparison.
#
# Usage: Rscript fog_deinvi6_figure.R [var_name ...]
#   var_name: one or more of fog (default: fog)

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
library(terra)
library(sf)
library(scatterpie)
library(ggspatial)

base_dir     <- "/private/groups/meyerlab/esteiger/selection/gemma"
env_file     <- "/private/groups/meyerlab/esteiger/env_var.txt"
field_file   <- "/private/groups/meyerlab/esteiger/field_data.txt"
scaffold_map <- file.path(base_dir, "scaffold_snp_map.txt")
dem_file     <- "/private/groups/meyerlab/esteiger/env/tranquillon_mountain_30m.dem"

box_fill <- "grey60"

base_theme <- theme_classic() +
  theme(legend.position  = "none",
        axis.title       = element_text(size = 22),
        axis.text        = element_text(size = 19),
        axis.line        = element_line(color = "black"),
        panel.border     = element_blank())

# Per-variable env-column name and axis label for the left panel. Both
# variables use the same "<var>_germ_yr" column naming convention in
# env_var.txt.
VAR_CONFIG <- list(
  fog = list(long_term_col = "fog", germ_yr_col = "fog_germ_yr",
             y_label = "Fog (days per month)")
)

args <- commandArgs(trailingOnly = TRUE)
var_names <- if (length(args) >= 1) args else names(VAR_CONFIG)

# ── Build the Tranquillon Mountain basemap once (shared across variables) ────
# CRS comes straight off the DEM's own header (UTM zone 10 / NAD27) rather
# than being hardcoded, then reprojected to lon/lat to line up with the
# Lat/Long columns in env_var.txt. Contour lines only (no hillshade raster).
dem    <- rast(dem_file)
dem_ll <- project(dem, "EPSG:4326")
contour_sf <- st_as_sf(as.contour(dem_ll, nlevels = 15))

# OHANA_order (stand-of-plants grouping) lives in field_data.txt, not
# env_var.txt — join it in by SampleID.
field_dt <- fread(field_file, select = c("SampleID", "OHANA_order"))
field_dt$SampleID <- as.character(field_dt$SampleID)

# Zoomed to the Strauss individuals themselves (range of their Lat/Long,
# +buffer) rather than the full DEM extent.
STRAUSS_BUFFER <- 0.003   # degrees
PIE_RADIUS     <- 0.0007  # degrees; one pie per OHANA_order stand
JITTER_RANGE   <- 0.0004  # degrees; tiny offset so co-located stands separate

make_map_panel <- function(stand_df, yr, map_xlim, map_ylim) {
  row <- stand_df[stand_df$Year == yr, ]
  row$r <- PIE_RADIUS
  ggplot() +
    geom_sf(data = contour_sf, color = "grey40", linewidth = 0.15,
            alpha = 0.7, inherit.aes = FALSE) +
    geom_scatterpie(data = row, aes(x = Long, y = Lat, r = r),
                     cols = c("standard", "inverted"), color = "black",
                     linewidth = 0.1, alpha = 0.85) +
    scale_fill_manual(values = c(standard = "white", inverted = "black")) +
    annotation_scale(location = "bl", width_hint = 0.25, style = "ticks",
                      line_width = 1, height = unit(0.15, "cm"),
                      text_cex = 0.7, pad_x = unit(0.3, "cm"),
                      pad_y = unit(0.3, "cm")) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE, datum = NA) +
    labs(title = yr) +
    theme_void() +
    theme(plot.title      = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.margin     = margin(2, 2, 2, 2))
}

# ── Identify the DEINVI_00006 peak (shared across variables) once ────────────
smap     <- fread(scaffold_map, select = c(1, 3), col.names = c("chr_num", "chr_name"))
smap     <- unique(smap)
chr6_num <- smap$chr_num[smap$chr_name == "DEINVI_00006"][1]
cat("DEINVI_00006 = internal chr", chr6_num, "\n")

# Select peak by position: the 69-SNP fog peak / 66-SNP BIO15 peak both
# start at ~48,86x,xxx (not the adjacent, smaller peak just below it).
TARGET_START <- 48864379L

for (var_name in var_names) {
  cfg <- VAR_CONFIG[[var_name]]
  if (is.null(cfg)) {
    cat("Skipping", var_name, "— no VAR_CONFIG entry (only fog/BIO15 supported).\n")
    next
  }
  cat("\n==========", var_name, "==========\n")

  out_dir       <- file.path(base_dir, var_name)
  long_term_col <- cfg$long_term_col
  germ_yr_col   <- cfg$germ_yr_col

  # ── Read outlier SNPs and assign peaks ───────────────────────────────────
  outlier_file <- file.path(out_dir, paste0("outliers_for_annot_40000_3_", var_name, ".txt"))
  outliers     <- fread(outlier_file, header = FALSE)
  snp_info <- data.frame(snp = as.character(outliers$V2),
                         chr = outliers$V1,
                         pos = outliers$V3,
                         stringsAsFactors = FALSE)
  snp_info <- snp_info[!duplicated(snp_info$snp), ]
  snp_info <- snp_info[order(snp_info$chr, snp_info$pos), ]

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

  chr6_snps   <- snp_info[snp_info$chr == chr6_num, ]
  pk_at_target <- chr6_snps$peak[which.min(abs(chr6_snps$pos - TARGET_START))]
  target_info  <- snp_info[snp_info$peak == pk_at_target, ]
  cat("Target peak:", pk_at_target, "on chr", chr6_num,
      "—", nrow(target_info), "SNPs, pos",
      min(target_info$pos), "–", max(target_info$pos), "\n")

  # ── Read dosages ──────────────────────────────────────────────────────────
  dosage_file <- file.path(out_dir, paste0("dosage_", var_name, ".txt"))
  ds          <- fread(dosage_file, header = TRUE)
  dosage      <- as.matrix(ds[, -1])
  rownames(dosage) <- as.character(ds$snp)
  dosage      <- t(dosage)

  # ── Read environmental data ───────────────────────────────────────────────
  env         <- fread(env_file, header = TRUE, sep = "\t")
  env$SampleID <- as.character(env$SampleID)
  env_match   <- match(rownames(dosage), env$SampleID)
  matched     <- !is.na(env_match)
  dosage      <- dosage[matched, , drop = FALSE]
  env_matched <- env[env_match[matched], ]

  long_term_env <- env_matched[[long_term_col]]
  germ_env      <- env_matched[[germ_yr_col]]
  sample_loc    <- as.character(env_matched$Location)
  sample_year   <- env_matched$Year

  # ── Polarize alleles (vectorized, controlling for coverage) ───────────────
  coverage <- as.numeric(env_matched$Cov)
  cov_mat  <- cbind(1, coverage)
  hat      <- cov_mat %*% solve(crossprod(cov_mat), t(cov_mat))
  env_resid <- as.numeric(long_term_env - hat %*% long_term_env)
  dos_resid <- dosage - hat %*% dosage
  betas     <- as.numeric(crossprod(dos_resid, env_resid) / colSums(dos_resid^2))
  flip      <- which(betas < 0)
  if (length(flip) > 0) dosage[, flip] <- 2 - dosage[, flip]
  cat("SNPs flipped:", length(flip), "of", ncol(dosage), "\n")

  # ── Subset to Strauss ───────────────────────────────────────────────────────
  strauss_idx    <- which(sample_loc == "Strauss")
  strauss_dosage <- dosage[strauss_idx, , drop = FALSE]

  df <- data.frame(
    germ_env = germ_env[strauss_idx],
    Year     = factor(sample_year[strauss_idx], levels = c("2020", "2021", "2022")),
    stringsAsFactors = FALSE
  )

  # ── Compute mean dosage for target peak ─────────────────────────────────────
  pk_cols   <- which(colnames(strauss_dosage) %in% target_info$snp)
  cat("Dosage columns found for target peak:", length(pk_cols), "\n")
  pk_score  <- rowMeans(strauss_dosage[, pk_cols, drop = FALSE], na.rm = TRUE)
  df$pk_score <- pk_score

  # ── Panel A: variable's germination-year values by year ─────────────────────
  p_env <- ggplot(df, aes(x = Year, y = germ_env)) +
    geom_boxplot(fill = box_fill, outlier.shape = NA, color = "black", width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.25, size = 1.6, color = "black") +
    labs(x = "Germination year",
         y = cfg$y_label) +
    base_theme

  # ── Panel B: mean favorable-allele dosage by year ────────────────────────────
  p_peak <- ggplot(df, aes(x = Year, y = pk_score)) +
    geom_boxplot(fill = box_fill, outlier.shape = NA, color = "black", width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.25, size = 1.6, color = "black") +
    labs(x = "Germination year",
         y = "Mean favorable allele dosage\nat putative inversion") +
    base_theme

  # ── Combine and save ──────────────────────────────────────────────────────────
  combined <- plot_grid(p_env, p_peak, ncol = 1, align = "v", axis = "lr")

  out_path <- file.path(out_dir, paste0(var_name, "_deinvi6_twopanel.jpg"))
  ggsave(out_path, combined, width = 5, height = 8, dpi = 200)
  cat("Saved:", out_path, "\n")

  # ── Alternate panel B: genotype-proportion stacked bar chart ─────────────────
  # pk_score isn't perfectly discrete (imputation/genotype-likelihood noise on
  # some individuals — see conversation), but SNPs in this peak are in tight
  # enough LD that it clusters right around 0/1/2 copies of the favorable
  # allele; round() to bin each individual into one of the 3 diploid genotype
  # classes at this locus.
  geno_df <- df %>%
    mutate(genotype = factor(round(pk_score), levels = c(0, 1, 2),
                              labels = c("Standard/Standard", "Heterozygous",
                                         "Inverted/Inverted"))) %>%
    count(Year, genotype, .drop = FALSE) %>%
    group_by(Year) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  p_geno <- ggplot(geno_df, aes(x = Year, y = prop, fill = genotype)) +
    geom_col(color = "black", width = 0.6, position = position_dodge2(padding = 0.1)) +
    scale_fill_manual(values = c("Standard/Standard" = "white",
                                  "Heterozygous"      = "grey50",
                                  "Inverted/Inverted" = "black")) +
    labs(x = "Germination year", y = "Genotype proportion") +
    base_theme

  combined_geno <- plot_grid(p_env, p_geno, ncol = 1, align = "v", axis = "lr")

  geno_out_path <- file.path(out_dir, paste0(var_name, "_deinvi6_genotype_by_year.jpg"))
  ggsave(geno_out_path, combined_geno, width = 5, height = 8, dpi = 200)
  cat("Saved:", geno_out_path, "\n")

  # ── 3-panel map: per-OHANA_order-stand haplotype frequency by year, on
  # the DEM basemap, zoomed to Strauss ──────────────────────────────────────
  # Each stand gets its own pie: "haplotype frequency" here is a 2-slice
  # allele frequency averaged across the stand's individuals: inverted =
  # mean(pk_score) / 2 (dosage is 0/1/2 copies of the favorable allele,
  # already polarized above per individual), standard = 1 - inverted.
  indiv_df <- df
  indiv_df$Lat         <- env_matched$Lat[strauss_idx]
  indiv_df$Long        <- env_matched$Long[strauss_idx]
  indiv_df$SampleID    <- env_matched$SampleID[strauss_idx]
  indiv_df$OHANA_order <- field_dt$OHANA_order[match(indiv_df$SampleID, field_dt$SampleID)]

  stand_df <- indiv_df %>%
    group_by(Year, OHANA_order) %>%
    summarise(Lat = mean(Lat), Long = mean(Long),
              inverted = mean(pk_score) / 2, n = n(), .groups = "drop") %>%
    mutate(standard = 1 - inverted)
  cat("\n--- DEINVI_00006 haplotype frequency by stand/year (", var_name, ") ---\n")
  print(stand_df)

  # Same small jitter per stand across all three years, so a stand doesn't
  # jump around panel-to-panel — just nudged enough to separate stands that
  # would otherwise plot on top of each other.
  set.seed(1)
  jitter_tbl <- data.frame(
    OHANA_order = unique(stand_df$OHANA_order),
    dLat  = runif(length(unique(stand_df$OHANA_order)), -JITTER_RANGE, JITTER_RANGE),
    dLong = runif(length(unique(stand_df$OHANA_order)), -JITTER_RANGE, JITTER_RANGE)
  )
  stand_df <- stand_df %>%
    left_join(jitter_tbl, by = "OHANA_order") %>%
    mutate(Lat = Lat + dLat, Long = Long + dLong)

  map_xlim <- range(stand_df$Long) + c(-1, 1) * STRAUSS_BUFFER
  map_ylim <- range(stand_df$Lat)  + c(-1, 1) * STRAUSS_BUFFER

  map_panels <- lapply(levels(df$Year), make_map_panel,
                        stand_df = stand_df, map_xlim = map_xlim, map_ylim = map_ylim)
  map_combined <- plot_grid(plotlist = map_panels, nrow = 1)

  map_out_path <- file.path(out_dir, paste0(var_name, "_deinvi6_map3panel.jpg"))
  ggsave(map_out_path, map_combined, width = 11.14, height = 5, dpi = 200)
  cat("Saved:", map_out_path, "\n")

  # ── Quick stats ───────────────────────────────────────────────────────────────
  cat("\n---", var_name, "germ-yr ANOVA by Year ---\n")
  print(summary(aov(germ_env ~ Year, data = df)))
  cat("\n--- DEINVI_00006 dosage ANOVA by Year (", var_name, ") ---\n")
  print(summary(aov(pk_score ~ Year, data = df)))
  cat("\n--- DEINVI_00006 dosage ~ germ_env regression (", var_name, ") ---\n")
  print(summary(lm(pk_score ~ germ_env, data = df)))

  cat("\n--- DEINVI_00006 genotype distribution by Year (", var_name, ") ---\n")
  genotype <- factor(round(df$pk_score), levels = c(0, 1, 2),
                      labels = c("Standard/Standard", "Heterozygous", "Inverted/Inverted"))
  geno_tab <- table(Year = df$Year, Genotype = genotype)
  print(geno_tab)
  cat("\nProportions:\n")
  print(round(prop.table(geno_tab, margin = 1), 3))
  cat("\nChi-square test (Year x genotype):\n")
  print(chisq.test(geno_tab))
  cat("\nFisher's exact test (Year x genotype, simulated p-value):\n")
  print(fisher.test(geno_tab, simulate.p.value = TRUE, B = 1e5))

  # ── Directed test: is dosage higher specifically in the highest-fog year? ───
  # One-sided (not just "do years differ") since the hypothesis has a
  # direction: the favorable allele should be more common in the year with
  # more germination-season fog. High-fog year determined from the data
  # itself (mean germ_env per year), not hardcoded.
  year_means <- df %>% group_by(Year) %>% summarise(mean_germ_env = mean(germ_env), .groups = "drop")
  cat("\n--- Mean", germ_yr_col, "by Year (", var_name, ") ---\n")
  print(year_means)
  high_year <- as.character(year_means$Year[which.max(year_means$mean_germ_env)])
  cat("Highest-fog germination year:", high_year, "\n")

  cat("\n--- One-sided t-test: dosage in", high_year, "> dosage in other years (", var_name, ") ---\n")
  print(t.test(df$pk_score[df$Year == high_year], df$pk_score[df$Year != high_year],
               alternative = "greater"))

  cat("\n--- One-sided pairwise t-tests: dosage in", high_year, "> dosage in each other year (", var_name, ") ---\n")
  other_years <- setdiff(levels(df$Year), high_year)
  for (yr in other_years) {
    cat("\n", high_year, "vs", yr, ":\n")
    print(t.test(df$pk_score[df$Year == high_year], df$pk_score[df$Year == yr],
                 alternative = "greater"))
  }

  cat("\n--- One-sided version of dosage ~ germ_env regression (", var_name, ") ---\n")
  reg <- lm(pk_score ~ germ_env, data = df)
  two_sided_p <- summary(reg)$coefficients["germ_env", "Pr(>|t|)"]
  beta <- coef(reg)["germ_env"]
  one_sided_p <- if (beta > 0) two_sided_p / 2 else 1 - two_sided_p / 2
  cat("beta =", signif(beta, 3), "; one-sided p (H1: beta > 0) =", signif(one_sided_p, 3), "\n")
}
