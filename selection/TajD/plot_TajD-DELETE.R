#!/usr/bin/env Rscript

# Usage: Rscript plot_TajD.R <outlier_file> [output_dir] [tag]
args <- commandArgs(trailingOnly = TRUE)
outlier_file <- if (length(args) >= 1) args[1] else stop("Usage: Rscript plot_TajD.R <outlier_file> [output_dir] [tag]")
outdir <- if (length(args) >= 2) args[2] else dirname(outlier_file)
tag   <- if (length(args) >= 3) args[3] else ""

if (getRversion() < "4.4.0") `%||%` <- function(x, y) if (is.null(x)) y else x

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")
library(data.table)
library(dplyr)
library(ggplot2)
library(ggridges)
library(stats)
library(viridis)
library(patchwork)

env_var  <- basename(outdir)
# file_var is used for output filenames; env_var is used for variable lookups
file_var <- if (nchar(tag)) paste0(env_var, "_", tag) else env_var

# Map variable name to env.txt column index
var_names <- c("PercH20", "CHILI", "BIO12", "BIO15", "BIO18", "fog",
               "fog_germ_yr", "BIO12_germ_yr", "BIO15_germ_yr", "BIO18_germ_yr",
               "BCM_30y", "BCM_germ_yr", "deficit", "severity", "drought_mn", "moisture")
env_col <- match(env_var, var_names)
if (is.na(env_col)) stop(paste("Unknown variable:", env_var))

# Read population assignments and environmental data
fam <- read.table("/private/groups/meyerlab/esteiger/selection/gemma/postParalog_MAF_4x_GLimp_clean.fam",
                  header = FALSE)
cluster <- fam$V1

env <- read.table("/private/groups/meyerlab/esteiger/selection/gemma/env.txt", header = FALSE)
env_vals <- env[[env_col]]

pop_env <- data.frame(cluster = cluster, env_val = env_vals)
pop_means <- pop_env %>%
    group_by(cluster) %>%
    summarise(mean_env = mean(env_val, na.rm = TRUE)) %>%
    arrange(mean_env)

message("Per-population mean ", env_var, ":")
for (r in seq_len(nrow(pop_means))) {
    message("  ", pop_means$cluster[r], ": ", round(pop_means$mean_env[r], 3))
}

# Fixed comparison: CDFW (large census N) vs Strauss (small census N)
pop_low  <- "Strauss"
pop_high <- "CDFW"
extreme_pops <- c(pop_low, pop_high)
col_low  <- viridis(10, option = "D")[1]
col_high <- viridis(10, option = "D")[10]
extreme_colors <- setNames(c(col_low, col_high), extreme_pops)

# SNP-to-scaffold mapping
snp_map <- fread("/private/groups/meyerlab/esteiger/selection/gemma/scaffold_snp_map.txt",
                 col.names = c("gemma_chr", "snp_id", "scaffold", "pos_on_scaffold"))

# Read outlier SNPs
outs_raw <- read.table(outlier_file, sep = " ", header = FALSE)
if (nrow(outs_raw) == 0) { message("Empty outlier file"); quit(save = "no") }
names(outs_raw)[1:2] <- c("gemma_chr", "snp_id")
outs <- merge(outs_raw[, c("gemma_chr", "snp_id")], snp_map,
              by = c("gemma_chr", "snp_id"), all.x = TRUE)
outs <- outs[!is.na(outs$scaffold), ]
if (nrow(outs) == 0) { message("No SNPs matched scaffold map"); quit(save = "no") }

# Chromosome lengths for genome-wide x-axis ordering
chroms <- fread("/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked.fai",
                select = 1:2, col.names = c("chr", "len"))
chroms <- chroms %>% arrange(desc(len))
chroms$start <- cumsum(c(0, head(chroms$len, -1)))

# Find all per-population TajD files
tajd_files <- Sys.glob("/private/groups/meyerlab/esteiger/selection/TajD/4x_*.Diana.out")
if (length(tajd_files) == 0) { message("No TajD files found"); quit(save = "no") }

# Map TajD filenames to .fam cluster names
tajd_to_cluster <- c(
    "strauss"    = "Strauss",
    "cdfw"       = "CDFW",
    "arguello"   = "VAFB_Arguello",
    "lionhead"   = "VAFB_Lion_Head",
    "conception" = "Point_Conception",
    "dangermond" = "Dangermond"
)

# Process all populations
all_outlier_d <- list()
all_genome_d <- list()

for (tf in tajd_files) {
    pop_name <- sub("^4x_", "", sub("\\.Diana\\.out$", "", basename(tf)))
    pop_name_cap <- tajd_to_cluster[pop_name]

    if (is.na(pop_name_cap) || !(pop_name_cap %in% pop_means$cluster)) {
        message("  Skipping TajD file ", basename(tf), ": no matching cluster")
        next
    }

    td <- fread(tf)
    td$chr <- as.character(td$chr)
    td <- merge(td, chroms, by = "chr", all.x = FALSE)
    td$WinCenter <- (as.numeric(td$w_start) + as.numeric(td$w_end)) / 2

    all_genome_d[[pop_name_cap]] <- data.frame(
        TajimasD = td$TajimasD, population = pop_name_cap, stringsAsFactors = FALSE)

    setkey(td, chr, WinCenter)
    td[, row_num := .I]
    outs_key <- data.table(chr = outs$scaffold, WinCenter = as.numeric(outs$pos_on_scaffold))
    setkey(outs_key, chr, WinCenter)
    matched_rows <- unique(na.omit(td[outs_key, row_num, roll = "nearest", nomatch = NA]))
    outlier_idx <- logical(nrow(td))
    if (length(matched_rows) > 0) outlier_idx[matched_rows] <- TRUE

    if (sum(outlier_idx) > 0) {
        all_outlier_d[[pop_name_cap]] <- data.frame(
            TajimasD = td$TajimasD[outlier_idx], population = pop_name_cap, stringsAsFactors = FALSE)
        message("  ", pop_name_cap, ": ", sum(outlier_idx), " outlier windows")
    }
}

if (length(all_outlier_d) == 0) { message("No outlier windows found"); quit(save = "no") }

outlier_df <- bind_rows(all_outlier_d)
genome_df <- bind_rows(all_genome_d)

# Merge mean env values
outlier_df <- merge(outlier_df, pop_means, by.x = "population", by.y = "cluster")
genome_df <- merge(genome_df, pop_means, by.x = "population", by.y = "cluster")

# Order populations by mean env value
pop_order <- pop_means$cluster[order(pop_means$mean_env)]
outlier_df$population <- factor(outlier_df$population, levels = pop_order)
genome_df$population <- factor(genome_df$population, levels = pop_order)

plot_theme <- theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

# Subsample genome background for plotting only (full data used for KS tests above)
genome_df_plot <- genome_df %>%
    group_by(population) %>%
    slice_sample(n = 50000) %>%
    ungroup()

# Build color scale: grey for genome, viridis for populations
pop_colors <- setNames(viridis(length(pop_order), option = "D"), pop_order)
all_colors <- c("Genome" = "grey70", pop_colors)

# ============================================================
# All-population plots: genome vs outlier paired per population
# ============================================================

# Build paired long-format: each population has genome-wide and outlier distributions
genome_df_plot$type <- "Genome"
outlier_df$type     <- "Outlier"
paired_data <- rbind(
    genome_df_plot[, c("TajimasD", "population", "type")],
    outlier_df[, c("TajimasD", "population", "type")]
)
paired_data$population <- factor(paired_data$population, levels = pop_order)
paired_data$type       <- factor(paired_data$type, levels = c("Genome", "Outlier"))

type_colors      <- c("Genome" = "grey70",  "Outlier" = "steelblue")
type_line_colors <- c("Genome" = "grey40",  "Outlier" = "steelblue4")

# --- Paired violin + boxplot per population ---
jpeg(file.path(outdir, paste0("TajD_violin_", file_var, ".jpg")),
     width = 10, height = 5, units = "in", res = 300)
print(
    ggplot(paired_data, aes(x = population, y = TajimasD, fill = type)) +
    geom_violin(alpha = 0.6, scale = "width", position = position_dodge(0.85)) +
    geom_boxplot(width = 0.12, position = position_dodge(0.85),
                 outlier.size = 0.3, alpha = 0.8) +
    scale_fill_manual(values = type_colors) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = NULL, y = "Tajima's D", fill = NULL) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
)
dev.off()

# --- Faceted density per population ---
jpeg(file.path(outdir, paste0("TajD_density_", file_var, ".jpg")),
     width = 10, height = 6, units = "in", res = 300)
print(
    ggplot(paired_data, aes(x = TajimasD, fill = type, color = type)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~population, scales = "free_y") +
    scale_fill_manual(values = type_colors) +
    scale_color_manual(values = type_line_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = "Tajima's D", y = "Density", fill = NULL, color = NULL) +
    plot_theme +
    theme(legend.position = "bottom")
)
dev.off()

# --- Faceted ECDF per population ---
jpeg(file.path(outdir, paste0("TajD_ecdf_", file_var, ".jpg")),
     width = 10, height = 6, units = "in", res = 300)
print(
    ggplot(paired_data, aes(x = TajimasD, color = type, linetype = type)) +
    stat_ecdf(linewidth = 1.0) +
    facet_wrap(~population) +
    scale_color_manual(values = c("Genome" = "grey40", "Outlier" = "steelblue")) +
    scale_linetype_manual(values = c("Genome" = "dashed", "Outlier" = "solid")) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    labs(x = "Tajima's D", y = "Cumulative proportion", color = NULL, linetype = NULL) +
    plot_theme +
    theme(legend.position = "bottom")
)
dev.off()

# KS tests per population: outlier vs that population's own genome background
for (pop in pop_order) {
    genome_d_pop  <- genome_df$TajimasD[genome_df$population == pop]
    outlier_d_pop <- outlier_df$TajimasD[outlier_df$population == pop]
    if (length(outlier_d_pop) >= 2) {
        ks <- ks.test(outlier_d_pop, genome_d_pop)
        message("  KS test (outlier vs own genome) ", pop, ": D=",
                round(ks$statistic, 3), ", p=", signif(ks$p.value, 3))
    }
}

# ============================================================
# Extreme-population two-panel figure
# ============================================================

env_extreme <- pop_env %>% filter(cluster %in% extreme_pops)
env_extreme$cluster <- factor(env_extreme$cluster, levels = extreme_pops)

wtest_env <- wilcox.test(env_val ~ cluster, data = env_extreme)
wp_env <- ifelse(wtest_env$p.value < 2.2e-16, "p < 2.2e-16",
                 paste0("p = ", signif(wtest_env$p.value, 3)))

pop_labels <- c("Strauss" = "Tranquillon Mtn", "CDFW" = "Pt Arguello")

panel_a <- ggplot(env_extreme, aes(x = cluster, y = env_val, fill = cluster)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1.5, width = 0.5) +
    scale_fill_manual(values = extreme_colors) +
    scale_x_discrete(labels = pop_labels) +
    annotate("text", x = 1.5, y = max(env_extreme$env_val, na.rm = TRUE),
             label = wp_env, size = 4.5, fontface = "italic") +
    labs(x = NULL, y = env_var) +
    plot_theme + theme(legend.position = "none",
                       axis.text.x = element_text(angle = 30, hjust = 1))

# Panel B: per-population genome vs outlier paired violin for extreme pops
extreme_paired <- paired_data[paired_data$population %in% extreme_pops, ]
extreme_paired$population <- factor(extreme_paired$population, levels = extreme_pops)

panel_b <- ggplot(extreme_paired, aes(x = population, y = TajimasD, fill = type)) +
    geom_violin(alpha = 0.6, scale = "width", position = position_dodge(0.85)) +
    geom_boxplot(width = 0.12, position = position_dodge(0.85),
                 outlier.size = 0.3, alpha = 0.8) +
    scale_fill_manual(values = type_colors) +
    scale_x_discrete(labels = pop_labels) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = NULL, y = "Tajima's D", fill = NULL) +
    plot_theme + theme(legend.position = "bottom",
                       axis.text.x = element_text(angle = 30, hjust = 1))

jpeg(file.path(outdir, paste0("TajD_extreme_", file_var, ".jpg")),
     width = 8, height = 4, units = "in", res = 300)
print(panel_a + panel_b + plot_layout(widths = c(1, 2)))
dev.off()

# --- Ridgeline version: 4 distributions (genome + outlier for each extreme pop) ---
ridge_data <- rbind(
    data.frame(
        TajimasD   = genome_df_plot$TajimasD[genome_df_plot$population %in% extreme_pops],
        population = as.character(genome_df_plot$population[genome_df_plot$population %in% extreme_pops]),
        type = "Genome"
    ),
    data.frame(
        TajimasD   = outlier_df$TajimasD[outlier_df$population %in% extreme_pops],
        population = as.character(outlier_df$population[outlier_df$population %in% extreme_pops]),
        type = "Outlier"
    )
)
ridge_levels <- c(
    paste0(pop_high, " (Outlier)"), paste0(pop_high, " (Genome)"),
    paste0(pop_low,  " (Outlier)"), paste0(pop_low,  " (Genome)")
)
ridge_data$group <- paste0(ridge_data$population, " (", ridge_data$type, ")")
ridge_data$group <- factor(ridge_data$group, levels = ridge_levels)
ridge_colors <- setNames(
    c(col_high, "grey70", col_low, "grey70"),
    ridge_levels
)

panel_b_ridge <- ggplot(ridge_data, aes(x = TajimasD, y = group, fill = group)) +
    geom_density_ridges(alpha = 0.6, scale = 0.9, rel_min_height = 0.01) +
    scale_fill_manual(values = ridge_colors) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = "Tajima's D", y = NULL) +
    plot_theme + theme(legend.position = "none")

jpeg(file.path(outdir, paste0("TajD_extreme_violin_", file_var, ".jpg")),
     width = 8, height = 4, units = "in", res = 300)
print(panel_a + panel_b_ridge + plot_layout(widths = c(1, 2)))
dev.off()

# Export per-population outlier TajD values
outlier_df$mean_env <- pop_means$mean_env[match(outlier_df$population, pop_means$cluster)]
write.table(outlier_df[, c("population", "TajimasD", "mean_env")],
            file.path(outdir, paste0("outlierswD_bypop_", file_var, ".txt")),
            quote = FALSE, row.names = FALSE, sep = "\t")

message("Done: TajD plots for ", file_var)
