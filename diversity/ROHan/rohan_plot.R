#!/usr/bin/env Rscript
# =============================================================================
# rohan_plot.R — Summarise ROHan output and plot ROH by population
#
# Closely mirrors roh_plot.R (bcftools roh equivalent); key differences:
#   - Input: rohan_regions.txt has a header row and no leading tag/quality cols
#   - all_samples scoped to >1x coverage individuals actually run by ROHan
#   - Output files prefixed rohan_ instead of roh_
#
# Input:  rohan_regions.txt  (produced by collect_rohan.sh)
#         /private/groups/meyerlab/esteiger/field_data.txt
#         /private/groups/meyerlab/esteiger/env_var.txt
#
# Output: rohan_by_cluster_{2,4,6}x.jpg
#         rohan_prop_long_by_cluster_{2,4,6}x.jpg
# =============================================================================

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")
library(ggplot2)

setwd("/private/groups/meyerlab/esteiger/diversity/ROH/ROHan")

# ==============================================================================
# 1. Load ROH regions
#    rohan_regions.txt columns (with header, produced by collect_rohan.sh):
#    Sample | Chrom | Begin | End | ROH_length_bp | Validated_sites
#    Source: .max.hmmrohl.gz files (highest posterior probability threshold)
# ==============================================================================
roh <- read.table("rohan_regions.txt", header = TRUE, sep = "\t")
roh$Sample <- as.character(roh$Sample)
# Standardise length column name to match roh_plot.R internals
roh$Length_bp <- roh$ROH_length_bp

cat(sprintf("ROH segments loaded: %d across %d samples\n",
            nrow(roh), length(unique(roh$Sample))))

# ==============================================================================
# 2. Per-individual summary
# ==============================================================================
roh_summary <- do.call(rbind, lapply(split(roh, roh$Sample), function(d) {
  data.frame(
    SampleID       = d$Sample[1],
    n_ROH          = nrow(d),
    total_ROH_Mb   = sum(d$Length_bp) / 1e6,
    mean_ROH_kb    = mean(d$Length_bp) / 1e3,
    median_ROH_kb  = median(d$Length_bp) / 1e3
  )
}))

# Individuals with zero ROH segments are absent from rohan_regions.txt —
# add them, scoped to the >1x samples that ROHan actually ran on
field <- read.table("/private/groups/meyerlab/esteiger/field_data.txt",
                    header = TRUE, sep = "\t")
field$SampleID <- as.character(field$SampleID)

all_samples <- data.frame(
  SampleID = field$SampleID[field$SampleID != "4531" &
                             !is.na(field$Cov) & field$Cov > 1]
)
roh_summary <- merge(all_samples, roh_summary, by = "SampleID", all.x = TRUE)
roh_summary[is.na(roh_summary$n_ROH), c("n_ROH", "total_ROH_Mb",
                                          "mean_ROH_kb", "median_ROH_kb")] <- 0

# ==============================================================================
# 3. Merge with metadata (include Cov for per-plot coverage filtering)
# ==============================================================================
roh_summary <- merge(roh_summary, field[, c("SampleID", "Clusters", "Cov")],
                     by = "SampleID", all.x = TRUE)

cat(sprintf("Samples in analysis: %d\n", nrow(roh_summary)))

# ==============================================================================
# 4. Merge fog (long-term average) from env_var.txt
# ==============================================================================
env <- read.table("/private/groups/meyerlab/esteiger/env_var.txt",
                  header = TRUE, sep = "\t")
env$SampleID <- as.character(env$SampleID)

fog_by_sample <- env[, c("SampleID", "fog")]
roh_summary <- merge(roh_summary, fog_by_sample, by = "SampleID", all.x = TRUE)

# Cluster-mean fog for ordering and colouring
fog_cluster <- aggregate(fog ~ Clusters, data = roh_summary, FUN = mean, na.rm = TRUE)
fog_cluster <- fog_cluster[order(fog_cluster$fog), ]

roh_summary$Clusters <- factor(roh_summary$Clusters,
                                levels = fog_cluster$Clusters)
roh_summary$mean_fog <- fog_cluster$fog[match(roh_summary$Clusters,
                                               fog_cluster$Clusters)]

cat("\nCluster order (increasing fog):\n")
print(fog_cluster)
cat("\nTotal ROH summary by cluster:\n")
print(aggregate(total_ROH_Mb ~ Clusters, data = roh_summary,
                FUN = function(x) round(c(median = median(x), mean = mean(x),
                                          IQR = IQR(x)), 2)))

# ==============================================================================
# 5. Plotting params (match inbreeding_plot.R / roh_plot.R style)
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

# ==============================================================================
# 6. Plot: total ROH length per individual — one plot per coverage threshold
#    viridis direction = -1  →  purple = high fog, yellow = low fog
# ==============================================================================
for (cov_thresh in c(2, 4, 6)) {
  dat <- roh_summary[!is.na(roh_summary$Cov) & roh_summary$Cov > cov_thresh, ]
  outfile <- sprintf("rohan_by_cluster_%dx.jpg", cov_thresh)
  cat(sprintf("Plotting %s (%d samples >%dx)...\n", outfile, nrow(dat), cov_thresh))
  jpeg(outfile, height = 1200, width = 1100)
  print(
    ggplot(dat, aes(x = Clusters, y = total_ROH_Mb)) +
      geom_boxplot(aes(fill = mean_fog), col = "black") +
      scale_x_discrete(labels = cluster_labels) +
      scale_fill_viridis_c(direction = -1, name = "Mean fog\n(hrs/day)",
                           option = "viridis") +
      labs(y = "Total ROH length (Mb)") +
      gvtp_theme +
      theme(axis.text        = element_text(vjust = 0, hjust = 1),
            axis.title.x     = element_blank(),
            legend.text      = element_text(size = 40),
            legend.title     = element_text(size = 40))
  )
  dev.off()
}

# ==============================================================================
# 7. ROH length stratification — proportion of total ROH in long segments
#    "Long" defined as >= 1 Mb, corresponding to common ancestry within
#    roughly the last ~50 generations (expected length = 100/(2g) Mb)
# ==============================================================================
LONG_THRESH_BP <- 1e6  # 1 Mb

roh_long <- do.call(rbind, lapply(split(roh, roh$Sample), function(d) {
  data.frame(
    SampleID     = d$Sample[1],
    total_ROH_bp = sum(d$Length_bp),
    long_ROH_bp  = sum(d$Length_bp[d$Length_bp >= LONG_THRESH_BP])
  )
}))

roh_long <- merge(all_samples, roh_long, by = "SampleID", all.x = TRUE)
roh_long[is.na(roh_long$total_ROH_bp), c("total_ROH_bp", "long_ROH_bp")] <- 0

# Individuals with no ROH at all get NA — excluded from plot
roh_long$prop_long_ROH <- ifelse(roh_long$total_ROH_bp > 0,
                                  roh_long$long_ROH_bp / roh_long$total_ROH_bp,
                                  NA)

roh_long <- merge(roh_long, field[, c("SampleID", "Clusters", "Cov")], by = "SampleID", all.x = TRUE)
roh_long <- merge(roh_long, fog_by_sample, by = "SampleID", all.x = TRUE)
roh_long$Clusters <- factor(roh_long$Clusters, levels = fog_cluster$Clusters)
roh_long$mean_fog <- fog_cluster$fog[match(roh_long$Clusters, fog_cluster$Clusters)]

for (cov_thresh in c(2, 4, 6)) {
  dat <- roh_long[!is.na(roh_long$Cov) & roh_long$Cov > cov_thresh, ]
  outfile <- sprintf("rohan_prop_long_by_cluster_%dx.jpg", cov_thresh)
  cat(sprintf("\nProportion long ROH summary by cluster (>%dx):\n", cov_thresh))
  print(aggregate(prop_long_ROH ~ Clusters, data = dat,
                  FUN = function(x) round(c(median = median(x, na.rm = TRUE),
                                            mean   = mean(x, na.rm = TRUE)), 3)))
  cat(sprintf("Plotting %s (%d samples >%dx)...\n", outfile, nrow(dat), cov_thresh))
  jpeg(outfile, height = 1200, width = 1200)
  print(
    ggplot(dat, aes(x = Clusters, y = prop_long_ROH)) +
      geom_boxplot(aes(fill = mean_fog), col = "black") +
      scale_x_discrete(labels = cluster_labels) +
      scale_fill_viridis_c(direction = -1, name = "Mean fog\n(hrs/day)",
                           option = "viridis") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(y = "Proportion of total ROH\nin long segments (>1 Mb)") +
      gvtp_theme +
      theme(axis.text        = element_text(vjust = 0, hjust = 1),
            axis.title.x     = element_blank(),
            axis.title.y     = element_text(size = 50),
            legend.text      = element_text(size = 40),
            legend.title     = element_text(size = 40),
            plot.margin      = unit(c(1, 1, 1, 4), "cm"))
  )
  dev.off()
  cat(sprintf("Plot written to %s\n", outfile))
}
