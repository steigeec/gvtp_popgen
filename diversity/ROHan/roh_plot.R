#!/usr/bin/env Rscript
# =============================================================================
# roh_plot.R — Summarise ROHan output and plot ROH by population
#
# Adapted from the original bcftools-roh version of this script (same name,
# formerly at diversity/ROH/roh_plot.R) to run on ROHan's output instead.
# Structure and plots mirror the original as closely as possible so the two
# are visually comparable; output files are all "rohan_"-prefixed so they
# never collide with the original bcftools-based jpgs still in diversity/ROH/.
#
# ROHan has completed for all 423 samples (see gvtp_popgen/diversity/ROHan/launch_rohan.sh).
# Unlike the original script (where a sample's absence from roh_regions.txt
# legitimately means "zero ROH found"), the completed-sample set is read from
# rohan_summary_stats.txt (one row per sample ROHan has finished) rather than
# from field_data.txt, and only THOSE samples get zero-imputed when they have
# no ROH segments. Sample counts are printed throughout for reference.
#
# Input:  rohan_regions.txt        (produced by collect_rohan.sh)
#         rohan_summary_stats.txt  (produced by collect_rohan.sh)
#         /private/groups/meyerlab/esteiger/diversity/ROH/roh_regions.txt   (old, bcftools roh)
#         /private/groups/meyerlab/esteiger/field_data.txt
#         /private/groups/meyerlab/esteiger/env_var.txt
#
# Output (each produced for BOTH metrics: {total_ROH, prop_long} — see
# Y_METRICS below; "prop_long" = proportion of the GENOME in long (>1Mb)
# ROH, an F_ROH-style fraction isolating recent inbreeding, not to be
# confused with "proportion of an individual's ROH burden that is long"):
#         rohan_{total_ROH,prop_long}_by_cluster_{2,4}x.jpg
#         rohan_{total_ROH,prop_long}_by_cluster_4x_lat.jpg
#         rohan_{total_ROH,prop_long}_vs_coverage.jpg    (coverage confound check)
#         rohan_{total_ROH,prop_long}_vs_{BCM_30y,BCM_germ_yr,fog_30y,fog_germ_yr,BIO18,BIO18_germ_yr}_{individuals,strauss}.jpg
#         rohan_total_ROH_strauss_by_year.jpg  (Tranquillon Mtn, 2020/2021/2022, single panel;
#                                                styled to match inbreeding_violin_by_year.jpg)
#         rohan_prop_long_strauss_by_year.jpg  (same, but 3-panel: warm-qtr precip / fog /
#                                                prop_long_ROH, all by year, same individuals,
#                                                skinny-tall panels via patchwork)
#         rohan_{total_ROH,prop_long}_bars_Strauss.jpg    (Tranquillon Mtn by stand/OHANA_order;
#                                                            styled to match inbreeding_bars_Strauss.jpg)
#         rohan_vs_bcftools_total_ROH_comparison.jpg      (old vs new, total ROH only)
#
# Showing both metrics everywhere is temporary while we're still deciding
# which one (or both) to use going forward.
# =============================================================================

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")
library(ggplot2)
library(ggbeeswarm)
library(patchwork)
library(dunn.test)

setwd("/private/groups/meyerlab/esteiger/diversity/ROH/ROHan")

OLD_ROH_REGIONS <- "/private/groups/meyerlab/esteiger/diversity/ROH/roh_regions.txt"

# ==============================================================================
# 1. Metadata
# ==============================================================================
field <- read.table("/private/groups/meyerlab/esteiger/field_data.txt",
                    header = TRUE, sep = "\t")
field$SampleID <- as.character(field$SampleID)

env <- read.table("/private/groups/meyerlab/esteiger/env_var.txt",
                  header = TRUE, sep = "\t")
env$SampleID <- as.character(env$SampleID)
fog_by_sample <- env[, c("SampleID", "fog")]

# The ROHan-eligible cohort (>3x coverage, prep_samples.sh) is 423 samples —
# NOT all 842 in field_data.txt, most of which are below the coverage cutoff
# and were never candidates for ROH analysis at all. This is the correct
# denominator for "how much of the run is done" progress reporting below.
n_eligible <- nrow(read.table("samples_gt3x.tsv", header = FALSE, sep = "\t"))

# Genome size (autosomal scaffolds >=1Mb, same set ROHan's --auto restricts
# to) — needed to express long-ROH burden as a proportion of the genome
# (F_ROH-style, using only long tracts), analogous to how total_ROH_Mb is
# reported as absolute Mb.
FAI <- "/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked.fai"
fai <- read.table(FAI, header = FALSE, sep = "\t",
                  col.names = c("scaffold", "length", "offset", "linebases", "linewidth"))
scaffolds_1mb <- read.table("scaffolds_1mb_names.txt", header = FALSE, sep = "\t")$V1
GENOME_SIZE_BP <- sum(fai$length[fai$scaffold %in% scaffolds_1mb])

# 1 Mb ~ 50 generations at ~1 cM/Mb (standard assumption without a linkage
# map) — the literature-standard cutoff for "recent inbreeding" (Ceballos et
# al. 2018; the wild/non-model-species convention in Kardos et al. and the
# Taylor et al. low-coverage ROH paper this pipeline's rohmu calibration
# already follows). NOT a smaller data-driven cutoff — see prior discussion.
LONG_THRESH_BP <- 1e6
LONG_THRESH_LABEL <- sprintf("%g kb", LONG_THRESH_BP / 1e3)

# ==============================================================================
# 2. Load ROHan (new) ROH regions + the set of samples ROHan has FINISHED
#    rohan_regions.txt columns (header, produced by collect_rohan.sh):
#    Sample | Chrom | Begin | End | ROH_length_bp | Validated_sites
#    Source: .mid.hmmrohl.gz files (HMM run with the point-estimate theta —
#    ROHan's actual best guess; .min/.max use the theta credible interval's
#    lower/upper bounds instead and are collected separately, see section 13)
# ==============================================================================
roh <- read.table("rohan_regions.txt", header = TRUE, sep = "\t")
roh$Sample <- as.character(roh$Sample)
roh$Length_bp <- roh$ROH_length_bp

# rohan_summary_stats.txt has exactly one row per sample ROHan has completed
# (whether or not it had any ROH) — this, not field_data.txt, is the correct
# "completed" set for zero-imputation.
rohan_done <- read.table("rohan_summary_stats.txt", header = TRUE, sep = "\t")
rohan_done$Sample <- as.character(rohan_done$Sample)
completed_samples <- data.frame(SampleID = rohan_done$Sample)

cat(sprintf("\n*** ROHan progress: %d / %d samples completed as of %s ***\n\n",
            nrow(completed_samples), n_eligible, Sys.Date()))
cat(sprintf("ROH segments loaded: %d across %d samples (of the completed set)\n",
            nrow(roh), length(unique(roh$Sample))))

# ==============================================================================
# 3. Per-individual summary (new / ROHan) — zero-impute only within the
#    completed set, NOT the full 423-sample cohort
# ==============================================================================
roh_summary <- do.call(rbind, lapply(split(roh, roh$Sample), function(d) {
  data.frame(
    SampleID       = d$Sample[1],
    n_ROH          = nrow(d),
    total_ROH_Mb   = sum(d$Length_bp) / 1e6,
    long_ROH_Mb    = sum(d$Length_bp[d$Length_bp >= LONG_THRESH_BP]) / 1e6,
    mean_ROH_kb    = mean(d$Length_bp) / 1e3,
    median_ROH_kb  = median(d$Length_bp) / 1e3
  )
}))

roh_summary <- merge(completed_samples, roh_summary, by = "SampleID", all.x = TRUE)
roh_summary[is.na(roh_summary$n_ROH), c("n_ROH", "total_ROH_Mb", "long_ROH_Mb",
                                          "mean_ROH_kb", "median_ROH_kb")] <- 0

# Proportion of the GENOME (not of total ROH) that's in long ROH — an
# F_ROH-style fraction computed using only long tracts, directly comparable
# in spirit to total_ROH_Mb but isolating the recent-inbreeding signal.
roh_summary$prop_genome_long_ROH <- roh_summary$long_ROH_Mb * 1e6 / GENOME_SIZE_BP

# ==============================================================================
# 4. Merge with metadata (Cov for per-plot filtering, Year for the Strauss plot,
#    OHANA_order for the by-stand plot matching inbreeding_bars_Strauss.jpg)
# ==============================================================================
roh_summary <- merge(roh_summary, field[, c("SampleID", "Clusters", "Cov", "Lat", "Year", "OHANA_order")],
                     by = "SampleID", all.x = TRUE)

cat(sprintf("Samples in analysis: %d\n", nrow(roh_summary)))

# ==============================================================================
# 5. Merge fog (long-term average); order clusters by mean fog
# ==============================================================================
roh_summary <- merge(roh_summary, fog_by_sample, by = "SampleID", all.x = TRUE)

fog_cluster <- aggregate(fog ~ Clusters, data = roh_summary, FUN = mean, na.rm = TRUE)
fog_cluster <- fog_cluster[order(fog_cluster$fog), ]

roh_summary$Clusters <- factor(roh_summary$Clusters, levels = fog_cluster$Clusters)
roh_summary$mean_fog <- fog_cluster$fog[match(roh_summary$Clusters, fog_cluster$Clusters)]

cat("\nCluster order (increasing fog):\n")
print(fog_cluster)
cat("\nTotal ROH summary by cluster (all completed samples):\n")
print(aggregate(total_ROH_Mb ~ Clusters, data = roh_summary,
                FUN = function(x) round(c(median = median(x), mean = mean(x),
                                          IQR = IQR(x), n = length(x)), 2)))

# ==============================================================================
# 6. Plotting params (match roh_plot.R / inbreeding_plot.R style)
# ==============================================================================
gvtp_theme <- theme_bw() +
  theme(
    text                 = element_text(size = 70),
    axis.text            = element_text(angle = 45, hjust = 1),
    plot.margin          = unit(c(1, 1, 1, 1), "cm"),
    legend.key.height    = unit(1, "cm"),
    legend.key.width     = unit(1.5, "cm"),
    legend.box.background = element_rect(fill = "white", color = "black"),
    legend.spacing.y     = unit(0.1, "cm"),
    panel.background     = element_blank(),
    panel.grid           = element_blank(),
    plot.background      = element_blank(),
    legend.background    = element_blank(),
    panel.border         = element_rect(colour = "black", fill = NA, linewidth = 2)
  )

cluster_labels <- c(
  "CDFW"             = "Gaviota",
  "Dangermond"       = "Dangermond",
  "Point_Conception" = "Pt. Conception",
  "Strauss"          = "Tranquillon Mtn",
  "VAFB_Arguello"    = "Pt. Arguello",
  "VAFB_Lion_Head"   = "Lion's Head"
)

popcolors <- c(
  "CDFW"             = "#e69f00",
  "CDFW_MR"          = "#e69f00",
  "Strauss"          = "#009e73",
  "Point_Conception" = "#0072b2",
  "VAFB_Lion_Head"   = "#f0e442",
  "VAFB_Arguello"    = "#d55e00",
  "Dangermond"       = "#56b4e9"
)

# Per-stand (OHANA_order) palette, copied verbatim from inbreeding_plot.R so
# the by-stand plot below is directly visually comparable to
# inbreeding_bars_Strauss.jpg.
clustercolors <- c(
  "72" = "#EE8866",
  "73" = "#99DDFF", "74" = "#99DDFF",
  "75" = "#EEDD88", "78" = "#EEDD88", "79" = "#EEDD88",
  "80" = "#EEDD88", "81" = "#EEDD88", "82" = "#EEDD88", "83" = "#EEDD88",
  "85" = "#DDDDDD", "86" = "#DDDDDD", "87" = "#DDDDDD", "89" = "#DDDDDD",
  "88" = "#77AADD", "90" = "#77AADD", "91" = "#77AADD", "392" = "#77AADD", "486" = "#77AADD"
)
OHANA_LEVELS <- c("72", "73", "74", "75", "78", "79", "80", "81", "82", "83",
                  "85", "86", "87", "89", "88", "90", "91", "392", "486")

# Both metrics are shown in every plot below, for now (total ROH: cumulative
# autozygosity across all time depths; proportion of genome in long ROH:
# isolates the recent-inbreeding signal specifically — see the >1Mb/50-
# generation note above). Once the analysis direction is settled, most
# plots below will likely drop to just one of these two.
Y_METRICS <- list(
  list(var = "total_ROH_Mb", lab = "Total ROH length (Mb)", tag = "total_ROH"),
  list(var = "prop_genome_long_ROH",
       lab = sprintf("Proportion of genome\nin long ROH (>%s)", LONG_THRESH_LABEL),
       tag = "prop_long")
)

# ==============================================================================
# 7. Plot: ROH burden per individual — one plot per coverage threshold,
#    per metric (total ROH and proportion of genome in long ROH)
# ==============================================================================
for (m in Y_METRICS) {
  for (cov_thresh in c(2, 4)) {
    dat <- roh_summary[!is.na(roh_summary$Cov) & roh_summary$Cov > cov_thresh, ]
    outfile <- sprintf("rohan_%s_by_cluster_%dx.jpg", m$tag, cov_thresh)
    cat(sprintf("Plotting %s (%d samples >%dx)...\n", outfile, nrow(dat), cov_thresh))
    jpeg(outfile, height = 1200, width = 1100)
    print(
      ggplot(dat, aes(x = Clusters, y = .data[[m$var]])) +
        geom_boxplot(aes(fill = mean_fog), col = "black") +
        scale_x_discrete(labels = cluster_labels) +
        scale_fill_viridis_c(direction = -1, name = "Mean fog\n(hrs/day)",
                             option = "viridis") +
        labs(y = m$lab,
             subtitle = sprintf("%d/%d ROHan samples completed", nrow(completed_samples), n_eligible)) +
        gvtp_theme +
        theme(axis.text        = element_text(vjust = 0, hjust = 1),
              axis.title.x     = element_blank(),
              legend.text      = element_text(size = 40),
              legend.title     = element_text(size = 40))
    )
    dev.off()
  }
}

# ==============================================================================
# 8. Coverage confound check — does ROH burden track sequencing depth rather
#    than biology? Two possible failure modes: low coverage could UNDER-call
#    ROH (less power to detect true homozygous stretches) or OVER-call it
#    (allele dropout at low depth creates false-homozygous sites that look
#    like ROH). Either way, if a metric just tracks coverage, it isn't
#    trustworthy as an inbreeding signal. Checked for both metrics.
#
#    Uses AnalysisCov = min(raw Cov, 6), NOT raw Cov — prep_samples.sh
#    downsamples anything above 6x down to exactly 6x, so that's the actual
#    depth ROHan saw. Plotting against raw Cov would scatter many different
#    original depths (all downsampled to the same 6x) across the x-axis for
#    no reason; AnalysisCov is what could actually drive a detection-power
#    confound.
# ==============================================================================
roh_summary$AnalysisCov <- pmin(roh_summary$Cov, 6)
dat_cov <- roh_summary[!is.na(roh_summary$AnalysisCov), ]
cat(sprintf("\n*** Coverage confound check, n=%d completed samples ***\n", nrow(dat_cov)))

for (m in Y_METRICS) {
  yvals <- dat_cov[[m$var]]
  if (nrow(dat_cov) >= 3 && sd(dat_cov$AnalysisCov) > 0 && sd(yvals) > 0) {
    r_cov <- cor(dat_cov$AnalysisCov, yvals, method = "pearson")
    cat(sprintf("Pearson correlation (analysis coverage vs %s): r = %.3f\n", m$var, r_cov))
  } else {
    r_cov <- NA
    cat(sprintf("Not enough variation / samples yet for a %s correlation coefficient.\n", m$var))
  }

  outfile <- sprintf("rohan_%s_vs_coverage.jpg", m$tag)
  jpeg(outfile, units = "in", height = 6, width = 7, res = 300)
  print(
    ggplot(dat_cov, aes(x = AnalysisCov, y = .data[[m$var]], fill = Clusters)) +
      geom_smooth(method = "lm", se = TRUE, color = "grey50",
                  linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
      geom_point(color = "black", size = 2.5, shape = 21, alpha = 0.8) +
      scale_fill_manual(values = popcolors, labels = cluster_labels, name = "Cluster") +
      scale_y_continuous(limits = c(0, NA)) +
      labs(x = "Analysis coverage (x; downsampled to 6x cap)",
           y = m$lab,
           title = "Coverage confound check",
           subtitle = sprintf("n=%d%s", nrow(dat_cov),
                              if (!is.na(r_cov)) sprintf("; r=%.2f", r_cov) else "")) +
      theme_bw(base_size = 14)
  )
  dev.off()
  cat(sprintf("Plot written to %s\n", outfile))
}

# ==============================================================================
# 9. ROH burden — alternative x-axis ordered by latitude (4x only), per metric
# ==============================================================================
lat_cluster <- aggregate(Lat ~ Clusters,
                         data = roh_summary[!is.na(roh_summary$Lat), ],
                         FUN = mean, na.rm = TRUE)
lat_cluster <- lat_cluster[order(lat_cluster$Lat, decreasing = TRUE), ]

roh_summary$Clusters_lat <- factor(roh_summary$Clusters, levels = lat_cluster$Clusters)

dat_lat <- roh_summary[!is.na(roh_summary$Cov) & roh_summary$Cov > 4, ]
if (nrow(dat_lat) >= 1) {
  for (m in Y_METRICS) {
    outfile <- sprintf("rohan_%s_by_cluster_4x_lat.jpg", m$tag)
    jpeg(outfile, height = 1200, width = 1200)
    print(
      ggplot(dat_lat, aes(x = Clusters_lat, y = .data[[m$var]])) +
        geom_boxplot(aes(fill = Clusters_lat), col = "black") +
        scale_x_discrete(labels = cluster_labels) +
        scale_fill_manual(values = popcolors, guide = "none") +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        labs(y = m$lab,
             subtitle = sprintf("%d/%d ROHan samples completed", nrow(completed_samples), n_eligible)) +
        gvtp_theme +
        theme(axis.text        = element_text(vjust = 0, hjust = 1),
              axis.title.x     = element_blank(),
              axis.title.y     = element_text(size = 50),
              plot.margin      = unit(c(1, 1, 1, 4), "cm"))
    )
    dev.off()
    cat(sprintf("Plot written to %s\n", outfile))
  }
  cat("\nCluster order by latitude (N to S):\n")
  print(lat_cluster)
} else {
  cat("\nSkipping latitude plot — no completed samples >4x yet.\n")
}

# ==============================================================================
# 10. BCM_30y, BCM_germ_yr, and BIO18 (warmest-quarter precip.) vs total ROH
# ==============================================================================
roh_bcm <- merge(roh_summary,
                 env[, c("SampleID", "BCM_30y", "BCM_germ_yr", "fog_germ_yr",
                         "BIO18", "BIO18_germ_yr")],
                 by = "SampleID", all.x = TRUE)
dat_bcm <- roh_bcm[!is.na(roh_bcm$Cov) & roh_bcm$Cov > 20, ]

cat(sprintf("\nEnvironmental scatter dataset: %d individuals (>20x)\n", nrow(dat_bcm)))

strauss_col <- popcolors["Strauss"]

make_env_scatter <- function(dat, xvar, xlab, yvar, ylab, outfile, subset_strauss = FALSE) {
  if (subset_strauss) dat <- dat[dat$Clusters == "Strauss", ]
  if (nrow(dat) < 3) {
    cat(sprintf("  Skipping %s — only %d completed sample(s) so far (need >=3).\n", outfile, nrow(dat)))
    return(invisible(NULL))
  }
  jpeg(outfile, units = "in", height = 5, width = 5, res = 300)
  p <- if (subset_strauss) {
    ggplot(dat, aes(x = .data[[xvar]], y = .data[[yvar]])) +
      geom_smooth(method = "lm", se = TRUE, color = "grey50",
                  linetype = "dashed", linewidth = 0.8) +
      geom_point(color = "black", fill = strauss_col, size = 2, shape = 21, alpha = 0.7)
  } else {
    ggplot(dat, aes(x = .data[[xvar]], y = .data[[yvar]], fill = Clusters)) +
      geom_smooth(method = "lm", se = TRUE, color = "grey50",
                  linetype = "dashed", linewidth = 0.8, aes(group = 1)) +
      geom_point(color = "black", size = 2, shape = 21, alpha = 0.7) +
      scale_fill_manual(values = popcolors, guide = "none")
  }
  print(p + labs(x = xlab, y = ylab,
                 subtitle = sprintf("n=%d", nrow(dat))) +
          theme_bw(base_size = 24))
  dev.off()
  cat(sprintf("  Written to %s (n=%d)\n", outfile, nrow(dat)))
}

ENV_VARS <- list(
  list(xvar = "BCM_30y",       xlab = "BCM (30y avg)",                              tag = "BCM_30y"),
  list(xvar = "BCM_germ_yr",   xlab = "BCM (germination year)",                     tag = "BCM_germ_yr"),
  list(xvar = "fog",           xlab = "Fog (30y avg, hrs/day)",                     tag = "fog_30y"),
  list(xvar = "fog_germ_yr",   xlab = "Fog (germination year, hrs/day)",            tag = "fog_germ_yr"),
  list(xvar = "BIO18",         xlab = "Warmest-quarter precip. (BIO18, 30y avg)",   tag = "BIO18"),
  list(xvar = "BIO18_germ_yr", xlab = "Warmest-quarter precip. (BIO18, germination year)", tag = "BIO18_germ_yr")
)

for (m in Y_METRICS) {
  for (ev in ENV_VARS) {
    make_env_scatter(dat_bcm, ev$xvar, ev$xlab, m$var, m$lab,
                     sprintf("rohan_%s_vs_%s_individuals.jpg", m$tag, ev$tag))
    make_env_scatter(dat_bcm, ev$xvar, ev$xlab, m$var, m$lab,
                     sprintf("rohan_%s_vs_%s_strauss.jpg", m$tag, ev$tag), subset_strauss = TRUE)
  }
}

# ==============================================================================
# 11. Old (bcftools roh) vs New (ROHan) comparison
#     Only meaningful for samples ROHan has actually completed.
#     Not expecting the same absolute scale: different underlying models
#     (bcftools roh from called genotypes vs ROHan directly from BAM
#     genotype likelihoods). rohan_regions.txt now uses ROHan's mid
#     (point-estimate) threshold, its best single guess rather than a bound.
#     What's worth checking here is whether the two broadly agree on WHICH
#     samples/clusters carry relatively more or less ROH, not the magnitude.
# ==============================================================================
old_roh <- read.table(OLD_ROH_REGIONS, header = FALSE, sep = "\t",
                      col.names = c("tag", "Sample", "Chrom", "Start", "End",
                                    "Length_bp", "N_markers", "Quality"))
old_roh$Sample <- gsub("_rmdup\\.bam$", "", old_roh$Sample)

old_summary <- do.call(rbind, lapply(split(old_roh, old_roh$Sample), function(d) {
  data.frame(SampleID = d$Sample[1], old_n_ROH = nrow(d), old_total_ROH_Mb = sum(d$Length_bp) / 1e6)
}))
old_summary <- merge(data.frame(SampleID = unique(field$SampleID)), old_summary, by = "SampleID", all.x = TRUE)
old_summary[is.na(old_summary$old_n_ROH), c("old_n_ROH", "old_total_ROH_Mb")] <- 0

new_summary <- roh_summary[, c("SampleID", "n_ROH", "total_ROH_Mb", "Clusters", "Cov")]
names(new_summary)[names(new_summary) == "n_ROH"] <- "new_n_ROH"
names(new_summary)[names(new_summary) == "total_ROH_Mb"] <- "new_total_ROH_Mb"

# Inner join, restricted to the completed-ROHan set — this is the whole point:
# comparing only where we actually have both numbers.
comparison <- merge(new_summary, old_summary, by = "SampleID")

cat(sprintf("\n*** Old (bcftools roh) vs New (ROHan) comparison, n=%d completed samples ***\n",
            nrow(comparison)))
print(comparison[order(-comparison$new_total_ROH_Mb),
                  c("SampleID", "Clusters", "Cov", "old_total_ROH_Mb", "new_total_ROH_Mb", "old_n_ROH", "new_n_ROH")])

if (nrow(comparison) >= 3 && sd(comparison$old_total_ROH_Mb) > 0 && sd(comparison$new_total_ROH_Mb) > 0) {
  r <- cor(comparison$old_total_ROH_Mb, comparison$new_total_ROH_Mb, method = "pearson")
  cat(sprintf("Pearson correlation (old vs new total ROH Mb): r = %.3f\n", r))
} else {
  r <- NA
  cat("Not enough variation / samples yet for a correlation coefficient.\n")
}

jpeg("rohan_vs_bcftools_total_ROH_comparison.jpg", units = "in", height = 6.5, width = 8.5, res = 300)
print(
  ggplot(comparison, aes(x = old_total_ROH_Mb, y = new_total_ROH_Mb, fill = Clusters)) +
    geom_point(color = "black", size = 3, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = popcolors, labels = cluster_labels, name = "Cluster") +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    labs(x = "Total ROH (Mb) — bcftools roh (old)",
         y = "Total ROH (Mb) — ROHan, mid threshold (new)",
         title = "Old vs new ROH calls",
         subtitle = sprintf("n=%d completed ROHan samples of %d total%s",
                            nrow(comparison), n_eligible,
                            if (!is.na(r)) sprintf("; r=%.2f", r) else "")) +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(size = 16), plot.subtitle = element_text(size = 12))
)
dev.off()
cat("Comparison plot written to rohan_vs_bcftools_total_ROH_comparison.jpg\n")

# ==============================================================================
# 12. Tranquillon Mtn (Strauss) ROH by year — 2020 / 2021 / 2022
#     ROHan (new) results only
#     Styled to match inbreeding_violin_by_year.jpg (inbreeding_plot.R) exactly:
#     plain gray80 boxplot (no per-year color, no jitter overlay), same
#     dimensions, same gvtp_theme + y-axis-title hjust tweak, no title/subtitle
#     clutter — so this is directly visually comparable to that plot.
#     (total_ROH only here — prop_long gets the 3-panel climate-context
#     version in section 12c below, per request.)
# ==============================================================================
strauss <- roh_summary[roh_summary$Clusters == "Strauss" & !is.na(roh_summary$Year), ]
strauss <- strauss[strauss$Year %in% c(2020, 2021, 2022), ]
strauss$Year <- factor(strauss$Year, levels = c(2020, 2021, 2022))

cat(sprintf("\n*** Tranquillon Mtn (Strauss) by year, n=%d completed samples ***\n", nrow(strauss)))
print(table(strauss$Year))

year_panel_theme <- gvtp_theme + theme(axis.title.y = element_text(hjust = 0.3))

if (nrow(strauss) >= 1) {
  outfile <- "rohan_total_ROH_strauss_by_year.jpg"
  jpeg(outfile, height = 1100, width = 1100)
  print(
    ggplot(data = strauss) +
      geom_boxplot(aes(x = Year, y = total_ROH_Mb), fill = "gray80", color = "black",
                   outlier.size = 5, linewidth = 2) +
      labs(x = "Year", y = "Total ROH length (Mb)") +
      year_panel_theme
  )
  dev.off()
  cat(sprintf("Plot written to %s\n", outfile))
} else {
  cat("Skipping Strauss-by-year plot — no completed Strauss samples with a 2020-2022 Year yet.\n")
}

# ==============================================================================
# 12c. Two-panel version for prop_long_ROH: BIO12 (water-year precip.) and
#      prop_long_ROH, both by year, same Strauss individuals, same
#      box-and-whisker style. Skinny/tall panels combined via patchwork,
#      following the multi-panel pattern in inbreeding_plot.R
#      (soilmoist_inbreeding_by_year_Strauss.jpg).
#
#      Uses BIO12_germ_yr (germination-year-specific value), NOT the 30-year
#      climate normal (BIO12) — the normal is fixed per location and
#      wouldn't show any real signal across sampling years; only the
#      germ-yr version actually varies year to year.
# ==============================================================================
strauss_clim <- merge(strauss, env[, c("SampleID", "BIO12_germ_yr")],
                      by = "SampleID", all.x = TRUE)

cat(sprintf("\n*** Strauss 2-panel climate-context plot, n=%d completed samples ***\n",
            nrow(strauss_clim)))

# --- Statistical test: prop_genome_long_ROH ~ Year (Strauss) -----------------
# Kruskal-Wallis (non-parametric, matches inbreeding_plot.R's kw_strauss_year
# approach for NGSF ~ Year) since ROH proportions aren't expected to be
# normally distributed. Dunn test with Bonferroni correction for pairwise
# year comparisons if the omnibus test suggests it's worth following up.
if (nrow(strauss_clim) >= 3 && length(unique(strauss_clim$Year)) >= 2) {
  cat("\nMedian prop_genome_long_ROH by year (Strauss):\n")
  print(tapply(strauss_clim$prop_genome_long_ROH, strauss_clim$Year, median))

  kw_prop_long <- kruskal.test(prop_genome_long_ROH ~ Year, data = strauss_clim)
  cat("\n=== Kruskal-Wallis: Strauss prop_genome_long_ROH ~ Year ===\n")
  print(kw_prop_long)

  cat("\n=== Dunn test (Bonferroni): Strauss prop_genome_long_ROH ~ Year ===\n")
  dunn_prop_long <- dunn.test(strauss_clim$prop_genome_long_ROH, strauss_clim$Year,
                              method = "bonferroni")

  # --- Drought-year contrast: 2021 (extreme drought, lowest BIO12_germ_yr)
  # vs. 2020+2022 combined (non-drought) — Wilcoxon rank-sum, matching the
  # binary drought-vs-non-drought contrast used for NGSF in inbreeding_plot.R.
  strauss_clim$drought_group <- ifelse(strauss_clim$Year == "2021",
                                       "drought_2021", "non_drought_2020_2022")
  cat("\n=== Wilcoxon rank sum: Strauss prop_genome_long_ROH, 2021 (drought) vs 2020+2022 (non-drought) ===\n")
  wilcox_prop_long_drought <- wilcox.test(prop_genome_long_ROH ~ drought_group,
                                          data = strauss_clim)
  print(wilcox_prop_long_drought)
} else {
  cat("\nSkipping Kruskal-Wallis/Dunn/Wilcoxon tests — not enough samples or year groups yet.\n")
}

if (nrow(strauss_clim) >= 1) {
  # CANVAS IS NOW FIXED at 4500x2400 (matching the very first sizing request)
  # and will NOT be rescaled again when these values change — scaling the
  # canvas in lockstep with font/linewidth (as earlier revisions did) makes
  # every version look identical at any given viewing size, since it's a
  # pure proportional scale-up. Growing these values against a fixed canvas
  # is what actually makes them bigger relative to the plot.
  # panel border linewidth 8, box/whisker linewidth 8, font size 160.
  big_panel_theme <- year_panel_theme +
    theme(text = element_text(size = 160),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 8))

  p_bcm <- ggplot(strauss_clim, aes(x = Year, y = BIO12_germ_yr)) +
    geom_boxplot(fill = "gray80", color = "black", outlier.size = 11, linewidth = 8) +
    labs(x = "Year", y = "Water year precip.") +
    big_panel_theme

  p_prop_long <- ggplot(strauss_clim, aes(x = Year, y = prop_genome_long_ROH)) +
    geom_boxplot(fill = "gray80", color = "black", outlier.size = 11, linewidth = 8) +
    labs(x = "Year", y = "Prop. ROH >1Mb") +
    big_panel_theme

  outfile <- "rohan_prop_long_strauss_by_year.jpg"
  jpeg(outfile, height = 2400, width = 4500)
  print(p_bcm + p_prop_long)
  dev.off()
  cat(sprintf("Plot written to %s\n", outfile))
} else {
  cat("Skipping Strauss 3-panel climate-context plot — no completed samples with a 2020-2022 Year yet.\n")
}

# ==============================================================================
# 12b. Tranquillon Mtn (Strauss) ROH by stand (OHANA_order)
#      Styled to match inbreeding_bars_Strauss.jpg exactly: boxplot per stand,
#      filled with the same clustercolors palette, x-axis labels hidden (too
#      many stands to label legibly), same gvtp_theme.
# ==============================================================================
strauss_stand <- roh_summary[roh_summary$Clusters == "Strauss" & !is.na(roh_summary$OHANA_order), ]
strauss_stand$OHANA_order <- factor(as.character(strauss_stand$OHANA_order), levels = OHANA_LEVELS)

cat(sprintf("\n*** Tranquillon Mtn (Strauss) by stand, n=%d completed samples ***\n", nrow(strauss_stand)))
print(table(strauss_stand$OHANA_order))

if (nrow(strauss_stand) >= 1) {
  for (m in Y_METRICS) {
    outfile <- sprintf("rohan_%s_bars_Strauss.jpg", m$tag)
    jpeg(outfile, height = 1000, width = 1200)
    print(
      ggplot() +
        geom_boxplot(data = strauss_stand,
                     aes(x = OHANA_order, y = .data[[m$var]], fill = OHANA_order),
                     colour = "black", linewidth = 2, outlier.size = 5) +
        scale_fill_manual(values = clustercolors, guide = "none") +
        labs(x = "Stand", y = m$lab) +
        gvtp_theme +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(hjust = 0))
    )
    dev.off()
    cat(sprintf("Plot written to %s\n", outfile))
  }
} else {
  cat("Skipping Strauss-by-stand plot — no completed Strauss samples with an OHANA_order yet.\n")
}

# ==============================================================================
# 13. Per-individual mid estimate with min/max spread
#     rohan_roh_burden_by_threshold.txt (from collect_rohan.sh) has, per
#     sample, total ROH Mb computed three ways: mid (point estimate, the bar),
#     max (theta upper bound -> conservative/fewer ROH -> a LOWER total), and
#     min (theta lower bound -> lenient/more ROH -> a HIGHER total).
#
#     IMPORTANT: the min-based (lenient) value is frequently extreme — in the
#     current data it reaches >500 Mb (out of a 1339 Mb genome) for 25 of 56
#     samples, i.e. it's calling nearly the WHOLE GENOME as ROH. This lines up
#     with ROHan's own FAQ, which says these interval bounds tend to run
#     inflated. Don't read the upper whisker in plot (a) as a tight
#     statistical bound — it's better read as "how extreme the estimate gets
#     under an implausibly low background-heterozygosity assumption." Plot
#     (b) below, using only the (much more reasonable) max-based lower bound,
#     is likely the more practically useful of the two.
# ==============================================================================
burden <- read.table("rohan_roh_burden_by_threshold.txt", header = TRUE, sep = "\t")
burden$Sample <- as.character(burden$Sample)
burden <- merge(burden, roh_summary[, c("SampleID", "Clusters")],
                by.x = "Sample", by.y = "SampleID", all.x = TRUE)
burden <- burden[order(burden$Clusters, -burden$Total_ROH_Mb_mid), ]
burden$Sample <- factor(burden$Sample, levels = unique(burden$Sample))

cat(sprintf("\n*** Per-individual mid/min/max ROH spread, n=%d completed samples ***\n", nrow(burden)))
cat(sprintf("Samples where the lenient (min) estimate exceeds 500 Mb (>37%% of the genome): %d of %d\n",
            sum(burden$Total_ROH_Mb_min > 500), nrow(burden)))

# --- (a) Full spread: mid bar, whisker from max (conservative) to min (lenient) —
#     pseudo-log y-axis, since values span ~4 orders of magnitude (0 to >1000 Mb)
n_facet_cols <- min(3, length(unique(burden$Clusters)))
jpeg("rohan_individual_spread_full.jpg",
     units = "in", height = 3 + 2.5 * ceiling(length(unique(burden$Clusters)) / n_facet_cols),
     width = 12, res = 300)
print(
  ggplot(burden, aes(x = Sample)) +
    geom_col(aes(y = Total_ROH_Mb_mid), fill = "grey70", color = "black", linewidth = 0.2) +
    geom_linerange(aes(ymin = Total_ROH_Mb_max, ymax = Total_ROH_Mb_min),
                   color = "firebrick", linewidth = 0.5, alpha = 0.7) +
    facet_wrap(~Clusters, scales = "free_x", ncol = n_facet_cols,
              labeller = as_labeller(cluster_labels)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                       breaks = c(0, 1, 10, 100, 1000)) +
    labs(x = NULL, y = "Total ROH (Mb), pseudo-log scale",
         title = "Per-individual ROH: mid estimate (bar) with min/max spread (red line)",
         subtitle = sprintf("n=%d completed samples. Red line = conservative (max) to lenient (min) range —\nlenient bound often approaches whole-genome size (see console output); read cautiously.", nrow(burden))) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          plot.subtitle = element_text(size = 10))
)
dev.off()
cat("Plot written to rohan_individual_spread_full.jpg\n")

# --- (b) Conservative-only range: mid bar, whisker down to max (conservative) —
#     linear scale; both values stay in a comparable, interpretable range
jpeg("rohan_individual_spread_conservative.jpg",
     units = "in", height = 3 + 2.5 * ceiling(length(unique(burden$Clusters)) / n_facet_cols),
     width = 12, res = 300)
print(
  ggplot(burden, aes(x = Sample)) +
    geom_col(aes(y = Total_ROH_Mb_mid), fill = "grey70", color = "black", linewidth = 0.2) +
    geom_linerange(aes(ymin = Total_ROH_Mb_max, ymax = Total_ROH_Mb_mid),
                   color = "steelblue4", linewidth = 0.8) +
    facet_wrap(~Clusters, scales = "free_x", ncol = n_facet_cols,
              labeller = as_labeller(cluster_labels)) +
    labs(x = NULL, y = "Total ROH (Mb)",
         title = "Per-individual ROH: mid estimate (bar) down to conservative (max) bound",
         subtitle = sprintf("n=%d completed samples. Blue line = mid down to the theta-upper-bound (conservative) estimate.", nrow(burden))) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          plot.subtitle = element_text(size = 10))
)
dev.off()
cat("Plot written to rohan_individual_spread_conservative.jpg\n")
