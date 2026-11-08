#!/usr/bin/env Rscript
# Compare environmental conditions across all localities for variables
# that have peak-filtered (40000_3) GWAS outliers.
# Output: env_compare_sites_40000_3.jpg

BASE_DIR <- "/private/groups/meyerlab/esteiger"

.libPaths(file.path(BASE_DIR, "software/R_libs"))
library(data.table)
library(dplyr)
library(ggplot2)

GEMMA_DIR <- file.path(BASE_DIR, "selection/gemma")
OUTDIR    <- file.path(BASE_DIR, "gvtp_popgen/selection")

# ── Population assignments ────────────────────────────────────────────────────
fam <- read.table(file.path(GEMMA_DIR, "postParalog_MAF_4x_GLimp_clean.fam"),
                  header = FALSE, stringsAsFactors = FALSE)

# 13 individuals with untrustworthy coordinates (Location "Strauss" but genetic cluster elsewhere)
EXCLUDE <- c("4548", "4549", "4550", "4551", "4552", "4553", "4554",
             "4555", "4556", "4557", "4559", "4561", "4562")
keep    <- !as.character(fam$V2) %in% EXCLUDE
message("Excluding ", sum(!keep), " individuals; ", sum(keep), " retained.")

cluster <- fam$V1[keep]

# ── Environmental data (rows = individuals, cols = variables) ─────────────────
env <- read.table(file.path(GEMMA_DIR, "env.txt"), header = FALSE)
env <- env[keep, , drop = FALSE]

var_names <- c("PercH20", "CHILI", "BIO12", "BIO15", "BIO18", "fog",
               "fog_germ_yr", "BIO12_germ_yr", "BIO15_germ_yr", "BIO18_germ_yr",
               "BCM_30y", "BCM_germ_yr", "deficit", "severity", "drought_mn", "moisture")

var_labels <- c(
    "PercH20"       = "Percent water",
    "CHILI"         = "CHILI index",
    "BIO12"         = "Annual precipitation (mm)",
    "BIO15"         = "Precipitation seasonality",
    "BIO18"         = "Precipitation, warmest quarter (mm)",
    "fog"           = "Fog (days / month)",
    "fog_germ_yr"   = "Fog, germination year (days / month)",
    "BIO12_germ_yr" = "Annual precipitation, germination year (mm)",
    "BIO15_germ_yr" = "Precipitation seasonality, germination year",
    "BIO18_germ_yr" = "Precipitation, warmest quarter, germination year (mm)",
    "BCM_30y"       = "BCM 30-year average",
    "BCM_germ_yr"   = "BCM, germination year",
    "deficit"       = "Climate water deficit (mm)",
    "severity"      = "Drought severity",
    "drought_mn"    = "Drought months",
    "moisture"      = "Soil moisture (mm)"
)

# ── Find variables with non-empty 40000_3 outlier files ──────────────────────
present_vars <- Filter(function(vn) {
    f <- file.path(GEMMA_DIR, vn, paste0("outliers_for_annot_40000_3_", vn, ".txt"))
    file.exists(f) && file.info(f)$size > 0
}, var_names)

message(length(present_vars), " variables with 40000_3 outliers: ",
        paste(present_vars, collapse = ", "))
if (length(present_vars) == 0) stop("No 40000_3 outlier files found.")

# ── Colour / label scheme ────────────────────────────────────────────────────
site_labels <- c(
    "CDFW"             = "CDFW",
    "Dangermond"       = "Dangermond",
    "Point_Conception" = "Pt. Conception",
    "Strauss"          = "Strauss",
    "VAFB_Arguello"    = "Pt. Arguello",
    "VAFB_Lion_Head"   = "Lion's Head"
)

pop_colors <- c(
    "CDFW"             = "#e69f00",
    "CDFW_MR"          = "#e69f00",
    "Strauss"          = "#009e73",
    "Point_Conception" = "#0072b2",
    "VAFB_Lion_Head"   = "#f0e442",
    "VAFB_Arguello"    = "#d55e00",
    "Dangermond"       = "#56b4e9"
)

# Population display order: group VAFB sites; otherwise alphabetical by internal name
pop_order <- c("VAFB_Lion_Head", "VAFB_Arguello", "Strauss",
               "Dangermond", "Point_Conception", "CDFW")
pop_order <- intersect(pop_order, unique(cluster))

# ── Build long-format data frame ──────────────────────────────────────────────
env_long <- bind_rows(lapply(present_vars, function(vn) {
    col_idx <- match(vn, var_names)
    data.frame(
        population = cluster,
        value      = env[[col_idx]],
        variable   = vn,
        stringsAsFactors = FALSE
    )
}))
env_long <- env_long[!is.na(env_long$value), ]
env_long$population <- factor(env_long$population, levels = pop_order)
env_long$variable   <- factor(env_long$variable,   levels = present_vars)

# ── Plot ──────────────────────────────────────────────────────────────────────
var_labeller <- as_labeller(
    function(x) ifelse(x %in% names(var_labels), var_labels[x], x))

plot_theme <- theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95"),
          strip.text       = element_text(face = "bold", size = 9),
          axis.text.x      = element_text(angle = 45, hjust = 1),
          legend.position  = "none")

p <- ggplot(env_long, aes(x = population, y = value, fill = population)) +
    geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 0.4, alpha = 0.85) +
    scale_fill_manual(values = pop_colors) +
    scale_x_discrete(labels = site_labels) +
    facet_wrap(~ variable, scales = "free_y", labeller = var_labeller) +
    labs(x = NULL, y = NULL) +
    plot_theme

n_cols <- min(3, length(present_vars))
n_rows <- ceiling(length(present_vars) / n_cols)

out_file <- file.path(OUTDIR, "env_compare_sites_40000_3.jpg")
jpeg(out_file, width = 4 * n_cols, height = 4 * n_rows, units = "in", res = 300)
print(p)
dev.off()
message("Written: ", out_file)
