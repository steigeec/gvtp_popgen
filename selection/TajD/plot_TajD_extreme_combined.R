#!/usr/bin/env Rscript
# Generates combined multi-panel TajD comparison plots (CDFW vs Strauss)
# across all environmental variables in a single figure.
# Usage: Rscript plot_TajD_extreme_combined.R [tag]
#   tag (optional): e.g. "40000_3" to use filtered-peak outlier files;
#                   omit or leave empty for full Bonferroni outlier files.
# Two output files: one violin version, one ridgeline version.

args <- commandArgs(trailingOnly = TRUE)
tag         <- if (length(args) >= 1 && nchar(args[1]) > 0) args[1] else ""
file_suffix <- if (nchar(tag)) paste0("_", tag) else ""

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")
library(data.table)
library(dplyr)
library(ggplot2)
library(ggridges)
library(patchwork)
library(grid)
library(lme4)
library(lmerTest)

# Polished site labels
site_labels <- c(
    "CDFW"             = "Gaviota",
    "Dangermond"       = "Dangermond",
    "Point_Conception" = "Pt. Conception",
    "Strauss"          = "Tranquillon Mtn",
    "VAFB_Arguello"    = "Pt. Arguello",
    "VAFB_Lion_Head"   = "Lion's Head"
)

# Bang Wong color palette — consistent across all population plots
pop_colors <- c(
    "CDFW"             = "#e69f00",
    "CDFW_MR"          = "#e69f00",
    "Strauss"          = "#009e73",
    "Point_Conception" = "#0072b2",
    "VAFB_Lion_Head"   = "#f0e442",
    "VAFB_Arguello"    = "#d55e00",
    "Dangermond"       = "#56b4e9"
)


TAJD_DIR   <- "/private/groups/meyerlab/esteiger/selection/TajD"
GEMMA_DIR  <- "/private/groups/meyerlab/esteiger/selection/gemma"
GENOME_FAI <- "/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked.fai"
OUTDIR     <- GEMMA_DIR
MIN_WINDOWS  <- 5   # minimum matched outlier windows required per population
EXCLUDE_VARS <- c("moisture")  # variables to exclude regardless of data

var_names <- c("PercH20", "CHILI", "BIO12", "BIO15", "BIO18", "fog",
               "fog_germ_yr", "BIO12_germ_yr", "BIO15_germ_yr", "BIO18_germ_yr",
               "BCM_30y", "BCM_germ_yr", "deficit", "severity", "drought_mn", "moisture")

var_labels <- c(
    "fog"      = "Fog (days / month)",
    "BIO18"    = "Precipitation-warmest quarter (mm)",
    "BIO15"    = "Precipitation seasonality",
    "moisture" = "Soil moisture (mm)",
    "deficit"  = "Climate water deficit (mm)"
)
panel_labeller <- as_labeller(function(x) ifelse(x %in% names(var_labels), var_labels[x], x))

tajd_to_cluster <- c(
    "strauss"    = "Strauss",
    "cdfw"       = "CDFW",
    "arguello"   = "VAFB_Arguello",
    "lionhead"   = "VAFB_Lion_Head",
    "conception" = "Point_Conception",
    "dangermond" = "Dangermond"
)

target_pops <- c("Strauss", "CDFW")
axis_labels  <- site_labels

# Putative inversion region
INV_SCAFFOLD <- "DEINVI_00006"
INV_START    <- 43e6
INV_END      <- 50e6

# SNP-to-scaffold mapping
snp_map <- fread(file.path(GEMMA_DIR, "scaffold_snp_map.txt"),
                 col.names = c("gemma_chr", "snp_id", "scaffold", "pos_on_scaffold"))

# Positions in tight LD (r2 >= 0.5, ngsLD, within 300kb) with each variable's
# outlier SNPs -- see ld_scan_outlier_snps.py. An outlier SNP is often just a
# marker for a haplotype under selection rather than the causal site itself,
# so windows nearest a position in tight LD with an outlier SNP are folded
# into the same outlier-window match set as windows nearest the SNP itself.
# Kept per-variable (var_name column) throughout, never pooled, to preserve
# the same fog/BIO15/BIO18-specific restriction used elsewhere in this script.
LD_LINKED_DIR <- file.path(GEMMA_DIR, "ld_outlier_windows", "linked_snp_positions")
ld_linked_files <- Sys.glob(file.path(LD_LINKED_DIR, "*.linked.tsv"))
ld_linked_files <- ld_linked_files[file.info(ld_linked_files)$size > 0]
ld_linked <- if (length(ld_linked_files) > 0) {
    rbindlist(lapply(ld_linked_files, fread,
                     col.names = c("var_name", "scaffold", "seed_pos", "pos", "r2")))
} else {
    data.table(var_name = character(), scaffold = character(),
              seed_pos = numeric(), pos = numeric(), r2 = numeric())
}
message("Loaded ", nrow(ld_linked), " LD-linked (variable, position) pairs from ",
       length(ld_linked_files), " scaffold(s)")

# Chromosome list for filtering TajD windows
chroms <- fread(GENOME_FAI, select = 1:2, col.names = c("chr", "len"))

# Load and key TajD tables for target pops only
tajd_files  <- Sys.glob(file.path(TAJD_DIR, "4x_*.Diana.out"))
tajd_by_pop <- list()
for (tf in tajd_files) {
    pop_name     <- sub("^4x_", "", sub("\\.Diana\\.out$", "", basename(tf)))
    pop_name_cap <- tajd_to_cluster[pop_name]
    if (is.na(pop_name_cap)) next
    td <- fread(tf)
    td$chr <- as.character(td$chr)
    td <- merge(td, chroms, by = "chr", all.x = FALSE)
    td$WinCenter <- (as.numeric(td$w_start) + as.numeric(td$w_end)) / 2
    setkey(td, chr, WinCenter)
    td[, row_num := .I]
    tajd_by_pop[[pop_name_cap]] <- td
    message("Loaded ", pop_name_cap, ": ", nrow(td), " windows")
}

plot_theme <- theme_bw(base_size = 10) +
    theme(panel.grid.minor  = element_blank(),
          strip.background  = element_rect(fill = "white"),
          strip.text        = element_text(face = "bold"))

if (FALSE) {  # other Tajima's D plots/stats — not needed for the Delta D figures below

# Collect outlier TajD windows per variable per population
outlier_rows     <- list()
inv_outlier_rows <- list()   # windows whose outlier SNPs fall in the putative inversion
for (vn in var_names) {
    if (nchar(tag)) {
        outlier_file <- file.path(GEMMA_DIR, vn,
                                  paste0("outliers_for_annot_", tag, "_", vn, ".txt"))
    } else {
        outlier_file <- file.path(GEMMA_DIR, vn,
                                  paste0("unfiltered_outliers_for_annot_", vn, ".txt"))
    }
    if (vn %in% EXCLUDE_VARS) next
    if (!file.exists(outlier_file) || file.info(outlier_file)$size == 0) next
    outs_raw <- tryCatch(read.table(outlier_file, sep = " ", header = FALSE),
                         error = function(e) NULL)
    if (is.null(outs_raw) || nrow(outs_raw) == 0) next
    names(outs_raw)[1:2] <- c("gemma_chr", "snp_id")
    outs <- merge(outs_raw[, c("gemma_chr", "snp_id")], snp_map,
                  by = c("gemma_chr", "snp_id"), all.x = TRUE)
    outs <- outs[!is.na(outs$scaffold), ]
    if (nrow(outs) == 0) next

    outs_key <- data.table(chr = outs$scaffold,
                           WinCenter = as.numeric(outs$pos_on_scaffold))
    setkey(outs_key, chr, WinCenter)

    inv_outs_key <- outs_key[chr == INV_SCAFFOLD &
                              WinCenter >= INV_START & WinCenter <= INV_END]
    if (nrow(inv_outs_key) > 0) setkey(inv_outs_key, chr, WinCenter)

    pop_counts   <- integer(length(target_pops))
    names(pop_counts) <- target_pops
    pop_rows     <- list()
    inv_pop_rows <- list()

    for (pop in target_pops) {
        td <- tajd_by_pop[[pop]]
        if (is.null(td)) next
        matched <- unique(na.omit(
            td[outs_key, row_num, roll = "nearest", nomatch = NA]))
        pop_counts[pop] <- length(matched)
        if (length(matched) > 0)
            pop_rows[[pop]] <- data.frame(
                TajimasD   = td$TajimasD[matched],
                Watterson  = td$Watterson[matched],
                population = pop,
                var_name   = vn,
                chr        = td$chr[matched],
                WinCenter  = td$WinCenter[matched],
                stringsAsFactors = FALSE)
        if (nrow(inv_outs_key) > 0) {
            inv_matched <- unique(na.omit(
                td[inv_outs_key, row_num, roll = "nearest", nomatch = NA]))
            if (length(inv_matched) > 0)
                inv_pop_rows[[pop]] <- data.frame(
                    TajimasD   = td$TajimasD[inv_matched],
                    population = pop,
                    var_name   = vn,
                    stringsAsFactors = FALSE)
        }
    }

    if (all(pop_counts >= MIN_WINDOWS)) {
        outlier_rows <- c(outlier_rows, pop_rows)
        if (length(inv_pop_rows) > 0)
            inv_outlier_rows <- c(inv_outlier_rows, inv_pop_rows)
        message("Kept ", vn, " (", paste(target_pops, pop_counts, sep = "=", collapse = ", "), ")")
    } else {
        message("Dropped ", vn, " (", paste(target_pops, pop_counts, sep = "=", collapse = ", "), ")")
    }
}

outlier_df     <- bind_rows(outlier_rows)
outlier_inv_df <- bind_rows(inv_outlier_rows)
present_vars   <- intersect(var_names, unique(outlier_df$var_name))
message(length(present_vars), " variables kept: ", paste(present_vars, collapse = ", "))

# Genome background: one sample per population (not replicated per variable)
genome_data <- bind_rows(lapply(target_pops, function(pop) {
    td <- tajd_by_pop[[pop]]
    data.frame(
        TajimasD   = sample(td$TajimasD, min(50000, nrow(td))),
        population = pop,
        track      = "Genome",
        stringsAsFactors = FALSE)
}))

# Outlier data: use var_name as the track identifier
outlier_plot <- outlier_df
outlier_plot$track <- outlier_plot$var_name

plot_data <- bind_rows(
    outlier_plot[, c("TajimasD", "population", "track")],
    genome_data
)

# Track order: Genome at top, env vars stacked below (bottom = first level)
track_order <- c(rev(present_vars), "Genome")
plot_data$population <- factor(plot_data$population, levels = c("Strauss", "CDFW"))
plot_data$track      <- factor(plot_data$track, levels = track_order)

# x_lims no longer used — ridge plots auto-scale to show full tails

track_colors <- setNames(rep("grey60", length(track_order)), track_order)

track_label_fn <- function(x) ifelse(x %in% names(var_labels), var_labels[x], x)
pop_label_fn   <- site_labels

plot_theme <- theme_bw(base_size = 10) +
    theme(panel.grid.minor  = element_blank(),
          strip.background  = element_rect(fill = "white"),
          strip.text        = element_text(face = "bold"))

panel_w        <- 3.5
panel_h_violin <- 4.0
panel_h_ridge  <- 3.5

# --- Version 1: Boxplot ---
p_violin <- ggplot(plot_data, aes(x = track, y = TajimasD, fill = population)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = pop_colors) +
    scale_x_discrete(labels = track_label_fn) +
    facet_wrap(~ population, ncol = 2, labeller = as_labeller(pop_label_fn)) +
    labs(x = NULL, y = "Tajima's D") +
    plot_theme +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1),
          legend.position = "none")

jpeg(file.path(OUTDIR, paste0("TajD_extreme_all_violin", file_suffix, ".jpg")),
     width = panel_w * 2, height = panel_h_violin,
     units = "in", res = 300)
print(p_violin)
dev.off()
message("Written: TajD_extreme_all_violin", file_suffix, ".jpg")

# --- Version 2: Ridgeline plots ---
p_ridge <- ggplot(plot_data, aes(x = TajimasD, y = track, fill = track)) +
    geom_density_ridges(alpha = 0.65, scale = 2.2, rel_min_height = 0.01) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = track_colors, labels = track_label_fn) +
    scale_y_discrete(labels = track_label_fn,
                     expand = expansion(add = c(0.1, 0.5))) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    facet_wrap(~ population, ncol = 2, labeller = as_labeller(pop_label_fn)) +
    labs(x = "Tajima's D", y = NULL) +
    plot_theme +
    theme(legend.position = "none")

jpeg(file.path(OUTDIR, paste0("TajD_extreme_all_ridge", file_suffix, ".jpg")),
     width = panel_w * 2, height = panel_h_ridge,
     units = "in", res = 300)
print(p_ridge)
dev.off()
message("Written: TajD_extreme_all_ridge", file_suffix, ".jpg")

# ============================================================
# Genome-wide Tajima's D comparison across all populations
# ============================================================
genome_all <- bind_rows(lapply(names(tajd_by_pop), function(pop) {
    data.frame(TajimasD = tajd_by_pop[[pop]]$TajimasD, population = pop,
               stringsAsFactors = FALSE)
}))

pop_medians <- genome_all %>%
    group_by(population) %>%
    summarise(med = median(TajimasD, na.rm = TRUE)) %>%
    arrange(med)
genome_all$population <- factor(genome_all$population, levels = pop_medians$population)

genome_all_plot <- bind_rows(lapply(split(genome_all, genome_all$population), function(df) {
    df[sample.int(nrow(df), min(50000L, nrow(df))), , drop = FALSE]
}))

pop_colors_all <- pop_colors

jpeg(file.path(OUTDIR, "TajD_allpops_genome.jpg"),
     width = 5, height = 5, units = "in", res = 300)
print(
    ggplot(genome_all_plot, aes(x = population, y = TajimasD, fill = population)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
    scale_fill_manual(values = pop_colors_all) +
    scale_x_discrete(labels = site_labels) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = NULL, y = "Tajima's D") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
)
dev.off()
message("Written: TajD_allpops_genome.jpg")

# ============================================================
# Per-window Watterson theta: genome-wide vs outlier-matched
# ============================================================

# Genome-wide Watterson for all populations (subsampled)
watterson_genome_all <- bind_rows(lapply(names(tajd_by_pop), function(pop) {
    td <- tajd_by_pop[[pop]]
    idx <- sample.int(nrow(td), min(50000L, nrow(td)))
    data.frame(Watterson = td$Watterson[idx], population = pop,
               type = "Genome", stringsAsFactors = FALSE)
}))

# Order populations by median genome-wide Watterson theta
w_pop_order <- watterson_genome_all %>%
    group_by(population) %>%
    summarise(med = median(Watterson, na.rm = TRUE)) %>%
    arrange(med) %>%
    pull(population)
watterson_genome_all$population <- factor(watterson_genome_all$population, levels = w_pop_order)

# Watterson theta: genome background vs per-variable outliers, faceted by population
watterson_target_genome <- watterson_genome_all[
    watterson_genome_all$population %in% target_pops, ]

ov_target <- outlier_df[outlier_df$population %in% target_pops, ]
watterson_by_var <- bind_rows(
    data.frame(Watterson  = watterson_target_genome$Watterson,
               population = watterson_target_genome$population,
               track      = "Genome", stringsAsFactors = FALSE),
    data.frame(Watterson  = ov_target$Watterson,
               population = ov_target$population,
               track      = ov_target$var_name, stringsAsFactors = FALSE)
)
w_track_order  <- c("Genome", present_vars)
w_track_labels <- c("Genome" = "Genome\nbackground",
                    setNames(
                        ifelse(present_vars %in% names(var_labels),
                               var_labels[present_vars], present_vars),
                        present_vars))
watterson_by_var$track      <- factor(watterson_by_var$track,      levels = w_track_order)
watterson_by_var$population <- factor(watterson_by_var$population, levels = target_pops)
watterson_by_var            <- watterson_by_var[!is.na(watterson_by_var$track), ]

p_watt <- ggplot(watterson_by_var, aes(x = track, y = Watterson, fill = population)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
    scale_fill_manual(values = pop_colors) +
    scale_x_discrete(labels = w_track_labels) +
    facet_wrap(~ population, ncol = 2, labeller = as_labeller(pop_label_fn)) +
    labs(x = NULL, y = "Watterson theta (per 100-site window)") +
    plot_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

jpeg(file.path(OUTDIR, paste0("TajD_watterson_outlier", file_suffix, ".jpg")),
     width = 7, height = 5, units = "in", res = 300)
print(p_watt)
dev.off()
message("Written: TajD_watterson_outlier", file_suffix, ".jpg")

# ============================================================
# Ridge plot with putative inversion outlier windows overlaid
# ============================================================
if (nrow(outlier_inv_df) > 0) {
    inv_vars <- intersect(present_vars, unique(outlier_inv_df$var_name))
    message(length(inv_vars), " variable(s) have inversion-region outliers: ",
            paste(inv_vars, collapse = ", "))

    # Build stacked track order: for each variable that has inversion outliers,
    # interleave an inversion-specific sub-track immediately below it.
    inv_track_suffx <- "_inv"
    stacked_order <- unlist(lapply(rev(present_vars), function(vn) {
        if (vn %in% inv_vars) c(paste0(vn, inv_track_suffx), vn) else vn
    }))
    stacked_order <- c(stacked_order, "Genome")

    inv_plot_data <- outlier_inv_df
    inv_plot_data$track      <- paste0(inv_plot_data$var_name, inv_track_suffx)
    inv_plot_data$population <- factor(inv_plot_data$population, levels = c("Strauss", "CDFW"))
    inv_plot_data <- inv_plot_data[inv_plot_data$var_name %in% inv_vars, ]

    stacked_data <- bind_rows(plot_data, inv_plot_data[, c("TajimasD", "population", "track")])
    stacked_data$track <- factor(stacked_data$track, levels = stacked_order)

    stacked_label_fn <- function(x) {
        is_inv <- grepl(inv_track_suffx, x, fixed = TRUE)
        base   <- sub(inv_track_suffx, "", x, fixed = TRUE)
        label  <- ifelse(base %in% names(var_labels), var_labels[base], base)
        ifelse(is_inv, paste0(label, " (inversion)"), label)
    }

    stacked_colors <- setNames(
        ifelse(grepl(inv_track_suffx, stacked_order, fixed = TRUE), "firebrick",
               track_colors[sub(inv_track_suffx, "", stacked_order, fixed = TRUE)]),
        stacked_order
    )
    stacked_colors["Genome"] <- track_colors["Genome"]

    p_ridge_inv <- ggplot(stacked_data, aes(x = TajimasD, y = track, fill = track)) +
        geom_density_ridges(alpha = 0.65, scale = 2.2, rel_min_height = 0.01) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        scale_fill_manual(values = stacked_colors, labels = stacked_label_fn) +
        scale_y_discrete(labels = stacked_label_fn,
                         expand = expansion(add = c(0.1, 0.5))) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        facet_wrap(~ population, ncol = 2, labeller = as_labeller(pop_label_fn)) +
        labs(x = "Tajima's D", y = NULL,
             caption = paste0("Red: outlier windows in putative inversion (",
                              INV_SCAFFOLD, " ", INV_START/1e6, "–", INV_END/1e6, " Mb)")) +
        plot_theme +
        theme(legend.position = "none",
              plot.caption = element_text(size = 7, color = "grey40"))

    jpeg(file.path(OUTDIR, paste0("TajD_extreme_all_ridge_inversion", file_suffix, ".jpg")),
         width = panel_w * 2, height = panel_h_ridge + 0.5 + 0.4 * length(inv_vars),
         units = "in", res = 300)
    print(p_ridge_inv)
    dev.off()
    message("Written: TajD_extreme_all_ridge_inversion", file_suffix, ".jpg")
} else {
    message("No inversion-region outlier windows found — skipping inversion overlay plot")
}

# ============================================================
# Statistical summary: median TajD + bootstrap CI + KS test
# Genome background vs per-variable outlier sets, per population
# ============================================================

set.seed(42)
N_BOOT <- 10000

boot_median_ci <- function(x, n_boot = N_BOOT) {
    x <- x[!is.na(x)]
    meds <- replicate(n_boot, median(sample(x, length(x), replace = TRUE)))
    c(n      = length(x),
      median = median(x),
      ci_lo  = unname(quantile(meds, 0.025)),
      ci_hi  = unname(quantile(meds, 0.975)))
}

stat_rows   <- list()
genome_meds <- setNames(numeric(length(target_pops)), target_pops)

for (pop in target_pops) {
    pop_label   <- if (pop %in% names(site_labels)) site_labels[[pop]] else pop
    genome_vals <- plot_data$TajimasD[
        plot_data$population == pop & plot_data$track == "Genome"]
    genome_vals <- genome_vals[!is.na(genome_vals)]

    bg <- boot_median_ci(genome_vals)
    genome_meds[pop] <- bg["median"]

    stat_rows[[paste0(pop, "__Genome")]] <- data.frame(
        population  = pop_label,
        group       = "Genome background",
        n           = bg["n"],
        median_TajD = bg["median"],
        ci_lo_95    = bg["ci_lo"],
        ci_hi_95    = bg["ci_hi"],
        delta_TajD  = 0,
        ks_D        = NA_real_,
        ks_p        = NA_real_,
        stringsAsFactors = FALSE)

    for (vn in present_vars) {
        out_vals <- plot_data$TajimasD[
            plot_data$population == pop & plot_data$track == vn]
        out_vals <- out_vals[!is.na(out_vals)]
        if (length(out_vals) < 2) next
        ob  <- boot_median_ci(out_vals)
        ks  <- ks.test(out_vals, genome_vals)
        grp <- if (vn %in% names(var_labels)) var_labels[[vn]] else vn
        stat_rows[[paste0(pop, "__", vn)]] <- data.frame(
            population  = pop_label,
            group       = grp,
            n           = ob["n"],
            median_TajD = ob["median"],
            ci_lo_95    = ob["ci_lo"],
            ci_hi_95    = ob["ci_hi"],
            delta_TajD  = ob["median"] - bg["median"],
            ks_D        = unname(ks$statistic),
            ks_p        = ks$p.value,
            stringsAsFactors = FALSE)
    }
}

stats_df <- bind_rows(stat_rows)
rownames(stats_df) <- NULL

# Pooled-outlier delta D per population (all env vars combined)
pooled_summary <- vapply(target_pops, function(pop) {
    pop_label   <- if (pop %in% names(site_labels)) site_labels[[pop]] else pop
    all_out     <- plot_data$TajimasD[
        plot_data$population == pop & plot_data$track != "Genome"]
    all_out     <- all_out[!is.na(all_out)]
    pool_med    <- median(all_out, na.rm = TRUE)
    delta       <- pool_med - genome_meds[pop]
    c(pop = pop_label, n_outlier_windows = length(all_out),
      genome_median = round(genome_meds[pop], 4),
      pooled_outlier_median = round(pool_med, 4),
      delta_TajD = round(delta, 4))
}, character(5))
pooled_df <- as.data.frame(t(pooled_summary), stringsAsFactors = FALSE)

out_stats <- file.path(OUTDIR, paste0("TajD_stats_summary", file_suffix, ".txt"))
write.table(stats_df, out_stats, sep = "\t", quote = FALSE, row.names = FALSE)
message("Written: TajD_stats_summary", file_suffix, ".txt")

report_lines <- c(
    "",
    "========== Tajima's D Statistical Summary ==========",
    paste0("Tag: ", if (nchar(file_suffix)) file_suffix else "(unfiltered)"),
    "",
    "--- Per-variable results ---",
    capture.output(print(stats_df, digits = 4, row.names = FALSE)),
    "",
    "--- Pooled-outlier delta D per population ---",
    capture.output(print(pooled_df, row.names = FALSE)),
    "  (delta_TajD = median across all outlier windows − genome background median)",
    "===================================================="
)
message(paste(report_lines, collapse = "\n"))

# ============================================================
# Inversion vs non-inversion outlier comparison (BIO15 and fog)
# ============================================================

inv_vars_focus <- intersect(c("BIO15", "fog"), present_vars)

if (length(inv_vars_focus) > 0 && nrow(outlier_df) > 0 &&
    all(c("chr", "WinCenter") %in% names(outlier_df))) {

    is_inv_win <- outlier_df$chr == INV_SCAFFOLD &
                  outlier_df$WinCenter >= INV_START &
                  outlier_df$WinCenter <= INV_END

    inv_cmp_rows <- list()
    for (pop in target_pops) {
        pop_label   <- if (pop %in% names(site_labels)) site_labels[[pop]] else pop
        genome_vals <- plot_data$TajimasD[
            plot_data$population == pop & plot_data$track == "Genome"]
        genome_vals <- genome_vals[!is.na(genome_vals)]

        for (vn in inv_vars_focus) {
            vn_mask <- outlier_df$population == pop & outlier_df$var_name == vn

            inv_vals  <- outlier_df$TajimasD[vn_mask &  is_inv_win]
            ninv_vals <- outlier_df$TajimasD[vn_mask & !is_inv_win]
            inv_vals  <- inv_vals[!is.na(inv_vals)]
            ninv_vals <- ninv_vals[!is.na(ninv_vals)]

            make_row <- function(vals, group_label) {
                if (length(vals) < 2) return(NULL)
                ob  <- boot_median_ci(vals)
                ks  <- ks.test(vals, genome_vals)
                data.frame(
                    population  = pop_label,
                    variable    = if (vn %in% names(var_labels)) var_labels[[vn]] else vn,
                    group       = group_label,
                    n           = ob["n"],
                    median_TajD = ob["median"],
                    ci_lo_95    = ob["ci_lo"],
                    ci_hi_95    = ob["ci_hi"],
                    delta_TajD  = ob["median"] - median(genome_vals),
                    ks_D        = unname(ks$statistic),
                    ks_p        = ks$p.value,
                    stringsAsFactors = FALSE)
            }

            inv_cmp_rows[[paste(pop, vn, "inv")]]  <- make_row(inv_vals,  "Inversion outliers")
            inv_cmp_rows[[paste(pop, vn, "ninv")]] <- make_row(ninv_vals, "Non-inversion outliers")
        }
    }

    inv_cmp_df <- bind_rows(inv_cmp_rows)
    rownames(inv_cmp_df) <- NULL

    out_inv <- file.path(OUTDIR, paste0("TajD_inversion_comparison", file_suffix, ".txt"))
    write.table(inv_cmp_df, out_inv, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Written: TajD_inversion_comparison", file_suffix, ".txt")

    inv_report <- c(
        "",
        "========== Inversion vs Non-inversion Outlier Comparison ==========",
        paste0("Tag: ", if (nchar(file_suffix)) file_suffix else "(unfiltered)"),
        paste0("Inversion region: ", INV_SCAFFOLD, " ", INV_START/1e6, "-", INV_END/1e6, " Mb"),
        "",
        capture.output(print(inv_cmp_df, digits = 4, row.names = FALSE)),
        "==================================================================="
    )
    message(paste(inv_report, collapse = "\n"))
} else {
    message("Skipping inversion comparison (BIO15/fog not in present_vars or chr/WinCenter missing)")
}

}  # end if (FALSE) — other Tajima's D plots/stats

# ============================================================
# Delta Tajima's D vs. inbreeding coefficient / fog variability
# One point per population per category (Background genome / Outliers /
# Outliers in inversion), pooled across all environmental variables.
# Unlike present_vars above, the outlier set for each (population, variable)
# here is kept if that population alone has >= MIN_WINDOWS matched windows
# (no joint requirement across populations).
# ============================================================
NGSF_FILE   <- "/private/groups/meyerlab/esteiger/diversity/ngsF/ngsF_results.txt"
ENVVAR_FILE <- "/private/groups/meyerlab/esteiger/env_var.txt"
BAMLIST_DIR <- "/private/groups/meyerlab/esteiger/mapping/merged"

# Bamlists used to build each population's SAF/TajD track (authoritative
# sample membership — matches tajd_by_pop exactly; Dangermond has no TajD
# track and is excluded)
pop_bamlists <- c(
    "CDFW"             = "bamlist_4x_cdfw.txt",
    "Strauss"          = "bamlist_4x_strauss.txt",
    "Point_Conception" = "bamlist_4x_conception.txt",
    "VAFB_Lion_Head"   = "bamlist_4x_lionhead.txt",
    "VAFB_Arguello"    = "bamlist_4x_arguello.txt"
)
allpop_names <- intersect(names(pop_bamlists), names(tajd_by_pop))

# "Outliers" is every matched outlier window (including ones inside the
# inversion); "Putative inversion outliers" is the subset of those windows
# that lie entirely inside the putative inversion.
allpop_outlier_rows     <- list()
allpop_inv_outlier_rows <- list()
for (vn in var_names) {
    if (vn %in% EXCLUDE_VARS) next
    if (nchar(tag)) {
        outlier_file <- file.path(GEMMA_DIR, vn,
                                  paste0("outliers_for_annot_", tag, "_", vn, ".txt"))
    } else {
        outlier_file <- file.path(GEMMA_DIR, vn,
                                  paste0("unfiltered_outliers_for_annot_", vn, ".txt"))
    }
    if (!file.exists(outlier_file) || file.info(outlier_file)$size == 0) next
    outs_raw <- tryCatch(read.table(outlier_file, sep = " ", header = FALSE),
                         error = function(e) NULL)
    if (is.null(outs_raw) || nrow(outs_raw) == 0) next
    names(outs_raw)[1:2] <- c("gemma_chr", "snp_id")
    outs <- merge(outs_raw[, c("gemma_chr", "snp_id")], snp_map,
                  by = c("gemma_chr", "snp_id"), all.x = TRUE)
    outs <- outs[!is.na(outs$scaffold), ]
    if (nrow(outs) == 0) next

    outs_key <- data.table(chr = outs$scaffold,
                           WinCenter = as.numeric(outs$pos_on_scaffold))

    # Fold in LD-linked positions (r2 >= 0.5) for this variable's outlier
    # SNPs, so windows in tight LD with an outlier SNP are matched the same
    # as windows nearest the SNP itself
    ld_vn <- ld_linked[var_name == vn, .(chr = scaffold, WinCenter = as.numeric(pos))]
    outs_key <- unique(rbind(outs_key, ld_vn))
    setkey(outs_key, chr, WinCenter)

    for (pop in allpop_names) {
        td <- tajd_by_pop[[pop]]
        if (is.null(td)) next
        matched <- unique(na.omit(
            td[outs_key, row_num, roll = "nearest", nomatch = NA]))
        if (length(matched) == 0) next

        if (length(matched) >= MIN_WINDOWS)
            allpop_outlier_rows[[paste(pop, vn)]] <- data.frame(
                TajimasD = td$TajimasD[matched], population = pop, var_name = vn,
                win_id = matched, stringsAsFactors = FALSE)

        w_chr   <- td$chr[matched]
        w_start <- as.numeric(td$w_start[matched])
        w_end   <- as.numeric(td$w_end[matched])
        entirely_inside <- w_chr == INV_SCAFFOLD & w_start >= INV_START & w_end <= INV_END

        inside_matched <- matched[entirely_inside]
        if (length(inside_matched) > 0)
            allpop_inv_outlier_rows[[paste(pop, vn)]] <- data.frame(
                TajimasD = td$TajimasD[inside_matched], population = pop, var_name = vn,
                win_id = inside_matched, stringsAsFactors = FALSE)
    }
}
allpop_outlier_df     <- bind_rows(allpop_outlier_rows)
allpop_inv_outlier_df <- bind_rows(allpop_inv_outlier_rows)

genome_med_all <- vapply(tajd_by_pop[allpop_names],
                         function(td) median(td$TajimasD, na.rm = TRUE), numeric(1))

# Build the population x category delta-D table. If restrict_var is given,
# only outlier windows for that single environmental variable are pooled
# (used for the fog-only version of these plots); otherwise all variables
# are pooled together.
build_pop_delta_df <- function(restrict_var = NULL) {
    out_df <- allpop_outlier_df
    inv_df <- allpop_inv_outlier_df
    if (!is.null(restrict_var)) {
        out_df <- out_df[out_df$var_name == restrict_var, ]
        inv_df <- inv_df[inv_df$var_name == restrict_var, ]
    }
    rows <- list()
    for (pop in allpop_names) {
        bg_med <- genome_med_all[[pop]]

        out_vals <- out_df$TajimasD[out_df$population == pop]
        out_vals <- out_vals[!is.na(out_vals)]
        if (length(out_vals) >= MIN_WINDOWS)
            rows[[paste0(pop, "_out")]] <- data.frame(
                population = pop, category = "Outliers",
                delta_TajD = median(out_vals) - bg_med, stringsAsFactors = FALSE)

        inv_vals <- inv_df$TajimasD[inv_df$population == pop]
        inv_vals <- inv_vals[!is.na(inv_vals)]
        if (length(inv_vals) >= MIN_WINDOWS)
            rows[[paste0(pop, "_inv")]] <- data.frame(
                population = pop, category = "Putative inversion outliers",
                delta_TajD = median(inv_vals) - bg_med, stringsAsFactors = FALSE)
    }
    bind_rows(rows)
}

pop_delta_df     <- build_pop_delta_df()
pop_delta_df_fog <- build_pop_delta_df(restrict_var = "fog")

# ============================================================
# KS tests: is each population's outlier-window Tajima's D distribution
# distinct from that population's own genome-wide background? Run for
# every population with a Tajima's D track, for both outlier categories,
# pooled across all environmental variables (same pooling as pop_delta_df).
# ============================================================
ks_rows <- list()
for (pop in allpop_names) {
    bg_vals <- tajd_by_pop[[pop]]$TajimasD
    bg_vals <- bg_vals[!is.na(bg_vals)]

    out_vals <- allpop_outlier_df$TajimasD[allpop_outlier_df$population == pop]
    out_vals <- out_vals[!is.na(out_vals)]
    if (length(out_vals) >= MIN_WINDOWS) {
        ks <- ks.test(out_vals, bg_vals)
        ks_rows[[paste0(pop, "_out")]] <- data.frame(
            population = pop, category = "Outliers", n = length(out_vals),
            ks_D = unname(ks$statistic), ks_p = ks$p.value, stringsAsFactors = FALSE)
    }

    inv_vals <- allpop_inv_outlier_df$TajimasD[allpop_inv_outlier_df$population == pop]
    inv_vals <- inv_vals[!is.na(inv_vals)]
    if (length(inv_vals) >= MIN_WINDOWS) {
        ks <- ks.test(inv_vals, bg_vals)
        ks_rows[[paste0(pop, "_inv")]] <- data.frame(
            population = pop, category = "Putative inversion outliers", n = length(inv_vals),
            ks_D = unname(ks$statistic), ks_p = ks$p.value, stringsAsFactors = FALSE)
    }
}
ks_df <- bind_rows(ks_rows)
ks_df$site <- site_labels[ks_df$population]
ks_df <- ks_df[, c("site", "population", "category", "n", "ks_D", "ks_p")]

out_ks <- file.path(OUTDIR, paste0("TajD_KS_tests", file_suffix, ".txt"))
write.table(ks_df, out_ks, sep = "\t", quote = FALSE, row.names = FALSE)
message("Written: ", basename(out_ks))
message(paste(capture.output(print(ks_df[, -2], row.names = FALSE, digits = 4)), collapse = "\n"))

# Mean inbreeding coefficient (ngsF F) per population, restricted to the
# same samples used to build that population's TajD/SAF track
bam_to_pop <- bind_rows(lapply(allpop_names, function(pop) {
    bams <- readLines(file.path(BAMLIST_DIR, pop_bamlists[[pop]]))
    data.frame(sample_id = sub("_rmdup\\.bam$", "", bams), population = pop,
               stringsAsFactors = FALSE)
}))

ngsf <- fread(NGSF_FILE)
ngsf$sample_id <- sub("_rmdup\\.bam$", "", ngsf$BAM)
ngsf_by_pop <- merge(ngsf, bam_to_pop, by = "sample_id") %>%
    group_by(population) %>%
    summarise(mean_NGSF = mean(NGSF, na.rm = TRUE), .groups = "drop")

# Mean fog days/month, and its variability over recent decades (CV and
# detrended SD — see env/extract_env.R), per population, same sample set as
# the TajD/SAF track
envvar <- fread(ENVVAR_FILE)
envvar$sample_id <- as.character(envvar$SampleID)
fog_by_pop <- merge(envvar, bam_to_pop, by = "sample_id") %>%
    group_by(population) %>%
    summarise(mean_fog          = mean(fog, na.rm = TRUE),
              mean_fog_CV       = mean(fog_CV, na.rm = TRUE),
              mean_fog_detrend  = mean(fog_detrend, na.rm = TRUE), .groups = "drop")

# Same, for precipitation seasonality (BIO15_CV / BIO15_detrend — see
# env/extract_env.R)
BIO15_by_pop <- merge(envvar, bam_to_pop, by = "sample_id") %>%
    group_by(population) %>%
    summarise(mean_BIO15          = mean(BIO15, na.rm = TRUE),
              mean_BIO15_CV       = mean(BIO15_CV, na.rm = TRUE),
              mean_BIO15_detrend  = mean(BIO15_detrend, na.rm = TRUE), .groups = "drop")

# Same, for the long-term Kling fog record (Kling_detrend — see
# env/extract_env.R; an alternative, longer-baseline measure of detrended
# fog variability than fog_detrend's 2000-2022 MODIS-based series)
Kling_by_pop <- merge(envvar, bam_to_pop, by = "sample_id") %>%
    group_by(population) %>%
    summarise(mean_Kling_detrend  = mean(Kling_detrend, na.rm = TRUE), .groups = "drop")

# Background genome sits at delta = 0 by definition (see y = 0 reference
# line drawn below) — only the two outlier categories are plotted as points
category_shapes <- c(
    "Outliers"                    = 21,  # circle
    "Putative inversion outliers" = 24   # triangle
)

delta_y_lab <- expression(Delta~"Tajima's D")

delta_scatter_theme <- plot_theme +
    theme(legend.position  = "bottom",
          axis.title       = element_text(size = 13),
          axis.text        = element_text(size = 11))

build_delta_plot <- function(df, x_col, x_lab,
                             point_size = 3, stroke = 0.6,
                             axis_title_size = 13, axis_text_size = 11,
                             legend_text_size = 11, legend_key_size = 1.2,
                             show_legend = TRUE,
                             axis_lines_only = FALSE, axis_line_size = 0.8,
                             hline_size = 0.5, plot_margin = margin(5.5, 5.5, 5.5, 5.5)) {
    legend_theme <- if (!show_legend) theme(legend.position = "none") else theme()

    axis_theme <- if (axis_lines_only) {
        theme(panel.border       = element_blank(),
              axis.line.x.bottom = element_line(color = "black", linewidth = axis_line_size),
              axis.line.y.left   = element_line(color = "black", linewidth = axis_line_size))
    } else {
        theme()
    }

    ggplot(df, aes(x = .data[[x_col]], y = delta_TajD, fill = population, shape = category)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = hline_size) +
        geom_point(size = point_size, color = "black", stroke = stroke, alpha = 0.9) +
        scale_fill_manual(values = pop_colors, guide = "none") +
        scale_shape_manual(values = category_shapes, name = NULL) +
        guides(shape = if (show_legend) guide_legend(override.aes = list(fill = "grey50")) else "none") +
        labs(x = x_lab, y = delta_y_lab) +
        delta_scatter_theme +
        legend_theme +
        axis_theme +
        theme(axis.title      = element_text(size = axis_title_size),
              axis.text       = element_text(size = axis_text_size),
              legend.text     = element_text(size = legend_text_size),
              legend.key.size = unit(legend_key_size, "lines"),
              plot.margin     = plot_margin)
}

save_delta_plot <- function(p, name_stub, width = 5.5, height = 5) {
    jpeg(file.path(OUTDIR, paste0(name_stub, file_suffix, ".jpg")),
         width = width, height = height, units = "in", res = 300)
    print(p)
    dev.off()
    message("Written: ", name_stub, file_suffix, ".jpg")
}

make_delta_df <- function(delta_df, env_df) {
    d <- merge(delta_df, env_df, by = "population")
    d$category <- factor(d$category, levels = names(category_shapes))
    d
}

# Pooled across all environmental variables
delta_ngsf_df        <- make_delta_df(pop_delta_df, ngsf_by_pop)
delta_fogmean_df     <- make_delta_df(pop_delta_df, fog_by_pop)
delta_fogcv_df       <- make_delta_df(pop_delta_df, fog_by_pop)
delta_fogdetrend_df  <- make_delta_df(pop_delta_df, fog_by_pop)

save_delta_plot(build_delta_plot(delta_ngsf_df, "mean_NGSF",
                                 "Inbreeding coefficient"),
                "TajD_delta_vs_ngsF")
save_delta_plot(build_delta_plot(delta_fogmean_df, "mean_fog", "Fog (days/month)"),
                "TajD_delta_vs_fogMean")
save_delta_plot(build_delta_plot(delta_fogcv_df, "mean_fog_CV", "CV of fog days/month (%)"),
                "TajD_delta_vs_fogCV")
save_delta_plot(build_delta_plot(delta_fogdetrend_df, "mean_fog_detrend",
                                 "Detrended SD of fog days/month"),
                "TajD_delta_vs_fogDetrend")

# Fog-only version — outliers restricted to the "fog" GEMMA scan
delta_ngsf_fog_df       <- make_delta_df(pop_delta_df_fog, ngsf_by_pop)
delta_fogmean_fog_df    <- make_delta_df(pop_delta_df_fog, fog_by_pop)
delta_fogcv_fog_df      <- make_delta_df(pop_delta_df_fog, fog_by_pop)
delta_fogdetrend_fog_df <- make_delta_df(pop_delta_df_fog, fog_by_pop)

save_delta_plot(build_delta_plot(delta_ngsf_fog_df, "mean_NGSF",
                                 "Inbreeding coefficient"),
                "TajD_delta_vs_ngsF_fog")
save_delta_plot(build_delta_plot(delta_fogmean_fog_df, "mean_fog", "Fog (days/month)"),
                "TajD_delta_vs_fogMean_fog")
save_delta_plot(build_delta_plot(delta_fogcv_fog_df, "mean_fog_CV", "CV of fog days/month (%)"),
                "TajD_delta_vs_fogCV_fog")
save_delta_plot(build_delta_plot(delta_fogdetrend_fog_df, "mean_fog_detrend",
                                 "Detrended SD of fog days/month"),
                "TajD_delta_vs_fogDetrend_fog")

# BIO15-only version — outliers restricted to the "BIO15" GEMMA scan
pop_delta_df_BIO15      <- build_pop_delta_df(restrict_var = "BIO15")
delta_BIO15_BIO15_df <- make_delta_df(pop_delta_df_BIO15, BIO15_by_pop)

save_delta_plot(build_delta_plot(delta_BIO15_BIO15_df, "mean_BIO15_detrend",
                                 "Detrended SD of\nprecip. seasonality",
                                 point_size = 8, stroke = 2,
                                 axis_title_size = 28, axis_text_size = 27, legend_text_size = 25,
                                 show_legend = FALSE,
                                 axis_lines_only = TRUE, axis_line_size = 0.7, hline_size = 1,
                                 plot_margin = margin(5.5, 12, 5.5, 5.5)),
                "TajD_delta_vs_BIO15_BIO15")

# Kling-fog version — same fog-scan outlier windows as TajD_delta_vs_fogDetrend_fog
# above (pop_delta_df_fog), but plotted against the longer-baseline Kling_detrend
# instead of the 2000-2022 MODIS-based fog_detrend. Added alongside the existing
# fog_detrend plot, not a replacement for it.
delta_klingdetrend_fog_df <- make_delta_df(pop_delta_df_fog, Kling_by_pop)

save_delta_plot(build_delta_plot(delta_klingdetrend_fog_df, "mean_Kling_detrend",
                                 "Detrended SD of fog (Kling, long-term)"),
                "TajD_delta_vs_KlingDetrend_fog")

# ============================================================
# Mixed-effects test: does Delta Tajima's D scale with an environmental
# metric, and does that relationship differ between "Outliers" and
# "Putative inversion outliers"?
#
# The population-level medians above collapse each population/category to
# a single point (n = 5 populations), which can't support a real
# significance test and is sensitive to any one population. This instead
# uses every individual outlier window as a data point
# (delta_TajD = window TajimasD - that population's own genome-wide
# median), with population as a random intercept to account for windows
# within a population not being independent draws (the environmental
# metric only varies at the population level, so naive per-window
# regression would be badly pseudoreplicated).
#
# Restricted to outlier windows from the "fog" GEA scan only (matching the
# *_fog descriptive plots above), not pooled across all environmental
# variables: this model tests whether Delta D in windows already implicated
# in fog-associated selection scales with fog variability (or inbreeding),
# so the windows going into it need to be fog-associated outliers, not GEA
# outliers for unrelated variables like precipitation seasonality or soil
# chemistry.
# ============================================================

dedupe_windows <- function(df) {
    df %>% distinct(population, win_id, .keep_all = TRUE) %>% select(-var_name, -win_id)
}

window_df <- bind_rows(
    dedupe_windows(allpop_outlier_df[allpop_outlier_df$var_name == "fog", ]) %>%
        mutate(category = "Outliers"),
    dedupe_windows(allpop_inv_outlier_df[allpop_inv_outlier_df$var_name == "fog", ]) %>%
        mutate(category = "Putative inversion outliers")
)
window_df$delta_TajD <- window_df$TajimasD - genome_med_all[window_df$population]
window_df$category   <- factor(window_df$category, levels = names(category_shapes))

# Bootstrap 95% CI of the median delta_TajD for each population x category,
# from the underlying window-level values (same approach used for TajD_stats
# above) — this is what actually gets plotted: an honest per-point CI rather
# than a fitted regression ribbon across only 5 populations.
set.seed(42)
boot_median_ci <- function(x, n_boot = 2000) {
    x <- x[!is.na(x)]
    meds <- replicate(n_boot, median(sample(x, length(x), replace = TRUE)))
    c(median = median(x), ci_lo = unname(quantile(meds, 0.025)),
      ci_hi  = unname(quantile(meds, 0.975)))
}
window_ci_df <- window_df %>%
    group_by(population, category) %>%
    summarise(ci = list(boot_median_ci(delta_TajD)), .groups = "drop") %>%
    mutate(median = vapply(ci, `[[`, numeric(1), "median"),
           ci_lo  = vapply(ci, `[[`, numeric(1), "ci_lo"),
           ci_hi  = vapply(ci, `[[`, numeric(1), "ci_hi")) %>%
    select(-ci) %>%
    as.data.frame()

fit_delta_mixed_model <- function(env_df, env_col) {
    d <- merge(window_df, env_df[, c("population", env_col)], by = "population")
    names(d)[names(d) == env_col] <- "env_x"
    # Random slope for category (not just a random intercept) is required
    # here: category varies *within* each population, so without a random
    # slope any population-specific idiosyncrasy in the category contrast
    # gets absorbed into the window-level residual instead of the
    # population-level random effect, which understates the true
    # uncertainty on the env_x:category interaction (the term that answers
    # "does the inversion category respond more steeply than other
    # outliers?") — that comparison is fundamentally limited by n = 5
    # populations, and the model should reflect that.
    lmer(delta_TajD ~ env_x * category + (1 + category | population), data = d, REML = TRUE)
}

# Satterthwaite-df estimate (via lmerTest::contest1D) of a linear
# combination of fixed effects — e.g. the slope within one category, or the
# difference in slopes between categories. Using the model's own
# Satterthwaite df (rather than a fixed z = 1.96 Wald interval) matters here
# because effects tied to env_x are only informed by 5 independent
# populations, so their effective df is close to 3, not the ~17,500 window
# observations — a plain Wald CI would understate the true uncertainty.
mixed_contrast <- function(m, ...) {
    b <- lme4::fixef(m)
    v <- setNames(rep(0, length(b)), names(b))
    args <- list(...)
    for (nm in names(args)) v[nm] <- args[[nm]]
    ct <- lmerTest::contest1D(m, v, confint = TRUE)
    data.frame(estimate = ct[["Estimate"]], se = ct[["Std. Error"]], df = ct[["df"]],
               ci_lo = ct[["lower"]], ci_hi = ct[["upper"]], p_value = ct[["Pr(>|t|)"]])
}

plot_mixed_model <- function(m, env_df, env_col, x_lab, name_stub,
                             point_size = 3, stroke = 0.6,
                             axis_title_size = 13, axis_text_size = 11,
                             legend_text_size = 11, legend_key_size = 1.2,
                             legend_inside = FALSE, show_legend = TRUE,
                             plot_width = 5.5, plot_height = 5,
                             axis_lines_only = FALSE, axis_line_size = 0.8,
                             hline_size = 0.5, plot_margin = margin(5.5, 5.5, 5.5, 5.5)) {
    inv_term <- grep("^env_x:category", names(lme4::fixef(m)), value = TRUE)

    slope_outliers  <- mixed_contrast(m, env_x = 1)
    slope_inversion <- do.call(mixed_contrast, c(list(m = m), setNames(list(1, 1), c("env_x", inv_term))))
    slope_diff      <- do.call(mixed_contrast, c(list(m = m), setNames(list(1), inv_term)))

    slope_df <- bind_rows(
        cbind(term = "Outliers slope",             slope_outliers),
        cbind(term = "Inversion outliers slope",    slope_inversion),
        cbind(term = "Difference (interaction)",    slope_diff)
    )
    slope_df$term <- factor(slope_df$term, levels = rev(slope_df$term))

    # Population medians with bootstrap 95% CI whiskers (from window_ci_df),
    # positioned on the continuous environmental-metric axis. No fitted
    # line/ribbon — with only 5 populations a regression band across them
    # implies more precision than the data support; the slope test itself
    # lives in the coefficient table written to the summary txt file.
    pts <- merge(window_ci_df, env_df[, c("population", env_col)], by = "population")
    names(pts)[names(pts) == env_col] <- "env_x"
    pts$category <- factor(pts$category, levels = names(category_shapes))

    legend_theme <- if (!show_legend) {
        theme(legend.position = "none")
    } else if (legend_inside) {
        theme(legend.position    = c(0.98, 0.02),
              legend.justification = c(1, 0),
              legend.direction   = "vertical",
              legend.background = element_rect(fill = "white", color = "grey70"),
              legend.margin      = margin(4, 6, 4, 6))
    } else {
        theme(legend.position = "bottom")
    }

    axis_theme <- if (axis_lines_only) {
        theme(panel.border       = element_blank(),
              axis.line.x.bottom = element_line(color = "black", linewidth = axis_line_size),
              axis.line.y.left   = element_line(color = "black", linewidth = axis_line_size))
    } else {
        theme()
    }

    p_main <- ggplot(pts, aes(x = env_x, y = median, fill = population, shape = category)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = hline_size) +
        geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), color = "black", linewidth = 0.5) +
        geom_point(size = point_size, color = "black", stroke = stroke, alpha = 0.9) +
        scale_fill_manual(values = pop_colors, guide = "none") +
        scale_shape_manual(values = category_shapes, name = NULL) +
        guides(shape = if (show_legend) guide_legend(override.aes = list(fill = "grey50")) else "none") +
        labs(x = x_lab, y = delta_y_lab) +
        delta_scatter_theme +
        legend_theme +
        axis_theme +
        theme(axis.title    = element_text(size = axis_title_size),
              axis.text     = element_text(size = axis_text_size),
              legend.text   = element_text(size = legend_text_size),
              legend.key.size = unit(legend_key_size, "lines"),
              plot.margin   = plot_margin)

    save_delta_plot(p_main, name_stub, width = plot_width, height = plot_height)

    slope_df
}

mixed_summary_lines <- character(0)
run_mixed <- function(env_df, env_col, x_lab, name_stub, ...) {
    m        <- fit_delta_mixed_model(env_df, env_col)
    slope_df <- plot_mixed_model(m, env_df, env_col, x_lab, name_stub, ...)
    mixed_summary_lines <<- c(mixed_summary_lines,
        paste0("\n===== ", name_stub, " (", env_col, ") ====="),
        capture.output(print(summary(m))),
        paste0("isSingular: ", isSingular(m)),
        "",
        "Slope contrasts (Satterthwaite df, 95% CI):",
        capture.output(print(slope_df, row.names = FALSE)))
}

run_mixed(ngsf_by_pop, "mean_NGSF",         "Inbreeding coefficient",         "TajD_mixedmodel_vs_ngsF")
run_mixed(fog_by_pop,  "mean_fog_detrend",  "Detrended SD of\nfog days/month", "TajD_mixedmodel_vs_fogDetrend",
          point_size = 8, stroke = 2,
          axis_title_size = 28, axis_text_size = 27, legend_text_size = 25,
          show_legend = FALSE,
          axis_lines_only = TRUE, axis_line_size = 0.7, hline_size = 1,
          plot_margin = margin(5.5, 12, 5.5, 5.5))
# We also tried this mixed model against the long-term Kling fog record's
# detrended SD (mean_Kling_detrend) instead of the Werner FLCC-derived
# fog_detrend above, as a longer-baseline sensitivity check. Kling_detrend
# ranked the 5 populations quite differently from fog_detrend (e.g. Pt.
# Conception moved from mid-pack to the highest-variability population, and
# Lion's Head/Arguello -- clearly distinct under fog_detrend -- were nearly
# tied), and its narrow range across populations (0.030-0.063) produced
# uninformative, extremely wide-CI slope estimates. We concluded the two
# metrics disagree too much to treat Kling_detrend as corroborating evidence
# and dropped it, using Werner FLCC (fog_detrend) as the sole basis for the
# fog-variability mixed model.

out_mm <- file.path(OUTDIR, paste0("TajD_mixedmodel_summary", file_suffix, ".txt"))
writeLines(mixed_summary_lines, out_mm)
message("Written: ", basename(out_mm))

# ============================================================
# Ridge plot: all populations, BIO18/BIO15/fog outliers + their inversion
# subset, colored by population, genome background always grey. Reuses
# allpop_outlier_df / allpop_inv_outlier_df from the mixed-model section
# above (already gated per-population by MIN_WINDOWS).
# ============================================================

ridge_pops <- c("VAFB_Lion_Head", "VAFB_Arguello", "Strauss",
                "Point_Conception", "CDFW")
ridge_pops_avail <- intersect(ridge_pops, names(tajd_by_pop))

genome_ridge <- bind_rows(lapply(ridge_pops_avail, function(pop) {
    td <- tajd_by_pop[[pop]]
    data.frame(TajimasD = sample(td$TajimasD, min(50000, nrow(td))),
               population = pop, track = "Genome", alpha_val = 0.8,
               stringsAsFactors = FALSE)
}))

build_ridge_plot <- function(ridge_vars, name_stub, height) {
    outlier_ridge_raw <- allpop_outlier_df[allpop_outlier_df$var_name %in% ridge_vars &
                                            allpop_outlier_df$population %in% ridge_pops_avail, ]
    outlier_ridge <- data.frame(TajimasD = outlier_ridge_raw$TajimasD,
                                population = outlier_ridge_raw$population,
                                track = outlier_ridge_raw$var_name,
                                alpha_val = 0.4, stringsAsFactors = FALSE)

    inv_ridge_raw <- allpop_inv_outlier_df[allpop_inv_outlier_df$var_name %in% ridge_vars &
                                            allpop_inv_outlier_df$population %in% ridge_pops_avail, ]
    inv_ridge <- data.frame(TajimasD = inv_ridge_raw$TajimasD,
                            population = inv_ridge_raw$population,
                            track = paste0(inv_ridge_raw$var_name, "_inv"),
                            alpha_val = 1, stringsAsFactors = FALSE)

    ridge_data <- bind_rows(genome_ridge, outlier_ridge, inv_ridge)

    ridge_track_order <- c(unlist(lapply(rev(ridge_vars), function(vn) c(paste0(vn, "_inv"), vn))),
                           "Genome")
    ridge_data$track      <- factor(ridge_data$track, levels = ridge_track_order)
    ridge_data$population <- factor(ridge_data$population, levels = ridge_pops)
    ridge_data$fill_group <- ifelse(ridge_data$track == "Genome", "Genome",
                                    as.character(ridge_data$population))
    ridge_fill_values <- c(pop_colors, "Genome" = "grey60")

    ridge_label_fn <- function(x) {
        is_inv <- grepl("_inv$", x)
        base   <- sub("_inv$", "", x)
        label  <- ifelse(base %in% names(var_labels), var_labels[base], base)
        label  <- ifelse(base == "Genome", "Genome background", label)
        ifelse(is_inv, paste0(label, " (inversion)"), label)
    }

    p_ridge <- ggplot(ridge_data,
                      aes(x = TajimasD, y = track, fill = fill_group, alpha = I(alpha_val))) +
        geom_density_ridges(scale = 2.2, rel_min_height = 0.01, color = "black", linewidth = 0.2) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        scale_fill_manual(values = ridge_fill_values) +
        scale_y_discrete(labels = ridge_label_fn, expand = expansion(add = c(0.1, 0.5))) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        facet_wrap(~ population, ncol = 3, labeller = as_labeller(site_labels)) +
        labs(x = "Tajima's D", y = NULL) +
        plot_theme +
        theme(legend.position = "none")

    jpeg(file.path(OUTDIR, paste0(name_stub, file_suffix, ".jpg")),
         width = 12, height = height, units = "in", res = 300)
    print(p_ridge)
    dev.off()
    message("Written: ", name_stub, file_suffix, ".jpg")
}

build_ridge_plot(c("BIO18", "BIO15", "fog"), "TajD_extreme_all_ridge_inversion", height = 7)
build_ridge_plot("fog", "TajD_extreme_all_ridge_inversion_fogOnly", height = 4)

message("Done.")
