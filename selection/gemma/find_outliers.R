#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
pheno_num <- args[1]
var_name  <- args[2]

dir.create(var_name, showWarnings=FALSE, recursive=TRUE)

report_path <- file.path(var_name, paste0(var_name, "_report.txt"))
# sink(report_path, split=TRUE)  # temporarily disabled — do not overwrite output files

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")
library(bit64)
library(data.table)
library(dplyr)
library(ggplot2)
library(snpStats)
library(GWAS.utils)
library(stats)
library(qqman)

froggy_theme <- theme(
  axis.title = element_text(size=65),
  axis.text.y = element_text(size=65, angle=45, hjust=1),
  axis.text.x = element_blank(),
  plot.margin = unit(c(0, 0, 0, 0), "cm"),
  legend.position = "none",
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank()
)

dic <- read.table("scaffold_snp_map.txt")
dic <- unique(dic[,c("V1", "V3")])

# Canonical scaffold order + cumulative genome offsets for the whole-genome
# Manhattan plots, derived from the reference genome's actual scaffold
# lengths — NOT from max(Position) within whatever random/outlier subset of
# SNPs is being plotted (the old approach). That approach recomputed order
# separately for each Manhattan variant from a randomly-downsampled 50% of
# non-outlier SNPs (no set.seed()), so two scaffolds of similar length could
# swap places between runs, or even between the two plots in the same run.
# Basing it on the fixed .fai lengths instead makes the left-to-right
# scaffold order identical every time, for every variable and every plot.
FAI <- "/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked.fai"
fai_df <- read.table(FAI, header = FALSE,
                     col.names = c("scaffold", "length", "offset", "linebases", "linewidth"),
                     stringsAsFactors = FALSE)
scaf2len <- setNames(fai_df$length, fai_df$scaffold)
chrom_order <- setNames(dic[, c("V1", "V3")], c("Chromosome", "scaffold"))
chrom_order$Chromosome <- as.character(chrom_order$Chromosome)
chrom_order$length <- scaf2len[chrom_order$scaffold]
chrom_order <- chrom_order[order(-chrom_order$length), ]
chrom_order$start <- c(0, head(cumsum(chrom_order$length), -1))

gemma <- read.table(paste0("output/pheno", pheno_num, "_lmm_results.assoc.txt"), header=T)
names(gemma)[1:3] <- c("Chromosome", "SNP", "Position")
print(paste("No rows in input:", nrow(gemma)))

# Calculate lambda
chisq_obs <- qchisq(gemma$p_lrt, df=1, lower.tail=FALSE)
lambda <- median(chisq_obs, na.rm=TRUE) / qchisq(0.5, df=1)
print(paste("Lambda:", round(lambda, 3)))

if (lambda > 1.05) {
  gemma$p_final <- pchisq(chisq_obs / lambda, df=1, lower.tail=FALSE)
  print("Applied genomic control")
  # jpeg(file.path(var_name, "QQplot-comparison.jpg"), width=1200, height=600)
  # par(mfrow=c(1,2))
  # qq(gemma$p_lrt, main=paste("Uncorrected, lambda =", round(lambda, 3)))
  # qq(gemma$p_final, main="GC Corrected")
  # dev.off()
} else {
  gemma$p_final <- gemma$p_lrt
  print("No genomic control needed")
  # jpeg(file.path(var_name, "QQplot.jpg"), width=600, height=600)
  # qq(gemma$p_lrt, main=paste("lambda =", round(lambda, 3)))
  # dev.off()
}

# Bonferroni threshold
thresh_bonf <- 0.05 / nrow(gemma)
bonf_hits <- sum(gemma$p_final < thresh_bonf, na.rm=TRUE)
print(paste("Bonferroni hits:", bonf_hits))

# FDR
fdr_results <- p.adjust(gemma$p_final, method = "fdr")
fdr_hits <- sum(fdr_results <= 0.05, na.rm=TRUE)
print(paste("FDR hits:", fdr_hits))

if (bonf_hits == 0 && fdr_hits == 0) {
  print("No outliers from Bonferroni or FDR — skipping Manhattan plot")
  # sink()
  quit(save="no", status=0)
}

# Tag outliers
gemma$out <- "no"
gemma$out[gemma$p_final < thresh_bonf] <- "out"
print(paste("Checking Bonferroni outliers:", print(table(gemma$out))))

# Downsample non-outliers
data <- gemma
data_no <- data[data$out=="no",]
data_no <- data_no[sample(nrow(data_no), nrow(data_no)/2), ]

pks <- as.data.table(data[data$out=="out",])
# Write-out temporarily disabled while iterating on the Manhattan plots below.
# write.table(pks, file.path(var_name, paste0("unfiltered_outliers_for_annot_", var_name, ".txt")),
#             quote=F, row.names=F, col.names=F)

# Sliding-window retention: keep SNPs with >= 2 neighbours within 40kb
if (nrow(pks) > 0) {
  win <- 40000; step <- 20000
  win_table <- list()
  for (ch in unique(pks$Chromosome)) {
    chrlen <- max(pks[Chromosome==ch, Position])
    starts <- seq(1, chrlen, by = step)
    ends <- pmin(starts + win - 1, chrlen)
    dtw <- data.table(Chromosome = ch, start = starts, end = ends)
    setkey(dtw, start, end)
    win_table[[as.character(ch)]] <- dtw
  }
  wins <- rbindlist(win_table)
  setkey(pks, Position, Position)
  wins2 <- wins[, .(Chromosome, wstart = start, wend = end)]
  setnames(wins2, c("Chromosome","wstart","wend"))
  setkey(wins2, wstart, wend)
  win_counts <- list()
  for (ch in unique(wins2$Chromosome)) {
    w <- wins2[Chromosome == ch]
    s <- pks[Chromosome == ch]
    if (nrow(s) == 0) next
    setnames(s, "Position", "pos")
    s[, pos2 := pos]
    setkey(s, pos, pos2)
    setkey(w, wstart, wend)
    ol <- foverlaps(w, s, by.x = c("wstart", "wend"), by.y = c("pos", "pos2"), type = "any", nomatch = 0L)
    cnt <- ol[, .N, by = .(Chromosome, wstart, wend)]
    setnames(cnt, c("Chromosome", "start", "end", "N"))
    win_counts[[ch]] <- cnt
    setnames(s, "pos", "Position")
    s[, pos2 := NULL]
  }
  win_counts_dt <- rbindlist(win_counts)
  win_kept <- win_counts_dt[N >= 3]
  print(paste("Windows with >=2 SNPs within 40kb:", nrow(win_kept)))

  setnames(win_kept, c("start", "end"), c("wstart", "wend"), skip_absent = TRUE)

  pks[, pos2 := Position]
  setkey(pks, Chromosome, Position, pos2)
  setkey(win_kept, Chromosome, wstart, wend)

  retained <- unique(foverlaps(pks, win_kept, by.x = c("Chromosome", "Position", "pos2"),
                        by.y = c("Chromosome", "wstart", "wend"), nomatch = 0L))
  retained[, c("wstart", "wend", "N", "pos2") := NULL]
  pks[, pos2 := NULL]

  # Write-out temporarily disabled while iterating on the Manhattan plots below.
  # write.table(unique(retained), file.path(var_name, paste0("outliers_for_annot_40000_3_", var_name, ".txt")),
  #             quote=F, row.names=F, col.names=F)
  print(paste("Retained SNPs (multi-SNP peaks):", nrow(unique(retained))))
} else {
  print("No Bonferroni outliers — skipping window analysis")
  retained <- pks
}

# Helper: add cumulative genomic position columns to a Manhattan data frame,
# using the fixed reference-genome scaffold order (chrom_order) so this is
# identical across calls/runs rather than derived from whatever data is
# passed in.
add_positions <- function(man, chrom_order) {
  man$Chromosome <- as.character(man$Chromosome)
  man <- merge(man, chrom_order[, c("Chromosome", "start")], by = "Chromosome")
  man$pos        <- man$Position + man$start
  man$Chromosome <- factor(man$Chromosome, levels = chrom_order$Chromosome)
  man <- man %>% arrange(pos)
  list(man = man, chroms = chrom_order)
}

# Helper: write plain + scaffold-labeled Manhattan plots
write_manhattan <- function(man, dic, thresh_bonf, var_name, tag = "") {
  REPS  <- round(length(levels(man$Chromosome)) / 2) + 1
  fstem <- if (nchar(tag)) paste0(tag, "_", var_name) else var_name

  jpeg(file.path(var_name, paste0("Manhattan_plot_", fstem, ".jpg")), width = 2000, height = 1000)
  print(ggplot(man) +
    xlab("Position along the genome") + ylab("-log10(p-value)") +
    geom_point(aes(x = pos, y = -log10(p_final), col = Chromosome, size = out, shape = out)) +
    scale_color_manual(values = rep(c("grey36", "black"), REPS)) +
    geom_hline(yintercept = -log10(thresh_bonf), col = "darkblue", lwd = 3) +
    scale_size_manual(values  = c("out" = 3, "no" = 2)) +
    scale_shape_manual(values = c("out" = 19, "no" = 19)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0.03, 0.03)) +
    froggy_theme)
  dev.off()

  man_lab <- merge(man, dic, by.x = "Chromosome", by.y = "V1")
  jpeg(file.path(var_name, paste0("Manhattan-label_", fstem, ".jpg")), width = 2000, height = 1000)
  print(ggplot(man_lab) +
    xlab("Position along the genome") + ylab("-log10(p-value)") +
    geom_jitter(aes(x = pos, y = -log10(p_final), col = Chromosome, size = out, shape = out),
                height = 0.04, stroke = 1.5) +
    scale_color_manual(values = rep(c("grey36", "black"), REPS)) +
    geom_hline(yintercept = -log10(thresh_bonf), col = "darkblue", lwd = 3) +
    scale_size_manual(values  = c("out" = 5, "no" = 2)) +
    scale_shape_manual(values = c("out" = 5, "no" = 19)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0.03, 0.03)) +
    geom_text(data = subset(man_lab, out == "out"),
              aes(x = pos, y = -log10(p_final) + 1, label = V3),
              col = "black", size = 4) +
    froggy_theme)
  dev.off()
}

if (TRUE) {

# ── Manhattan 1: all Bonferroni outliers ──────────────────────────────────────
res_all <- add_positions(rbind(data_no, as.data.frame(pks)), chrom_order)
# TEMPORARILY DISABLED — regenerating only the per-scaffold 40000_3 plots below
# write_manhattan(res_all$man, dic, thresh_bonf, var_name)

# ── Manhattan 2: filtered peaks only (40000_3 set) ────────────────────────────
if (nrow(retained) > 0) {
  retained$out <- "out"
  res_ret <- add_positions(rbind(data_no, as.data.frame(retained)), chrom_order)
  # TEMPORARILY DISABLED — regenerating only the per-scaffold 40000_3 plots below
  # write_manhattan(res_ret$man, dic, thresh_bonf, var_name, tag = "40000_3")

  # ── Per-scaffold Manhattan plots for 40000_3 peaks (scaffolds > 40 Mb only) ──
  long_scaffolds <- fai_df$scaffold[fai_df$length > 40e6]

  # dic has V1 = gemma_chr, V3 = scaffold name (one scaffold per gemma chromosome)
  chr2scaf <- setNames(dic$V3, dic$V1)

  # Scaffold name for each retained outlier
  retained_df <- as.data.frame(retained)
  retained_df$scaffold <- chr2scaf[as.character(retained_df$Chromosome)]

  # Cluster retained SNPs into peak intervals (40kb merge gap, same logic as
  # select_peak_labels.R) so each peak can get a background band on its
  # scaffold's plot, padded a bit wider than the peak itself.
  PEAK_MERGE_GAP <- 40000
  retained_df <- retained_df[order(retained_df$scaffold, retained_df$Position), ]
  peak_id <- integer(nrow(retained_df)); cur <- 0; prev_scaf <- NA; prev_pos <- -Inf
  for (i in seq_len(nrow(retained_df))) {
    if (!identical(retained_df$scaffold[i], prev_scaf) || (retained_df$Position[i] - prev_pos) > PEAK_MERGE_GAP) cur <- cur + 1
    peak_id[i] <- cur
    prev_scaf <- retained_df$scaffold[i]; prev_pos <- retained_df$Position[i]
  }
  retained_df$peak_id <- peak_id
  peaks_by_scaffold <- retained_df %>%
    group_by(peak_id, scaffold) %>%
    summarise(Peak_Start = min(Position), Peak_End = max(Position), .groups = "drop")

  # Only process scaffolds with retained outliers that are also > 40 Mb
  target_scaffolds <- sort(intersect(
    unique(retained_df$scaffold[!is.na(retained_df$scaffold)]),
    long_scaffolds
  ))
  print(paste(length(target_scaffolds), "long scaffold(s) with 40000_3 outliers:",
              paste(target_scaffolds, collapse = ", ")))

  if (length(target_scaffolds) > 0) {
    # Add scaffold name to all GEMMA SNPs (data has Chromosome, SNP, Position, p_final, out)
    data$scaffold <- chr2scaf[as.character(data$Chromosome)]

    for (scaf in target_scaffolds) {
      # Points = non-outlier background SNPs (full resolution, not the
      # genome-wide-downsampled data_no) plus the 40000_3 retained peak SNPs
      # — this drops only the scattered Bonferroni hits that did NOT survive
      # the 40kb/>=3-neighbour window filter, while keeping the full
      # non-outlier background for context.
      non_outlier_scaf <- data[!is.na(data$scaffold) & data$scaffold == scaf & data$out == "no",
                                c("Position", "p_final")]
      retained_scaf <- retained_df[!is.na(retained_df$scaffold) & retained_df$scaffold == scaf,
                                    c("Position", "p_final")]
      scaf_snps <- rbind(non_outlier_scaf, retained_scaf)
      if (nrow(scaf_snps) == 0) next

      scaf_len <- fai_df$length[fai_df$scaffold == scaf]
      x_breaks <- seq(0, scaf_len, by = 10e6)
      x_labels <- paste0(x_breaks / 1e6, " Mb")

      # Background bands for this scaffold's peaks — wider than the peak
      # (and wider than the plotted points themselves) so it's clearly
      # visible, full plot height. A percentage-of-peak-width pad alone is
      # invisible/thinner than the points at whole-scaffold zoom (peaks are
      # often <100kb on scaffolds tens of Mb long), so floor it at 1.5% of
      # this scaffold's length per side instead of a fixed bp value.
      scaf_peaks <- peaks_by_scaffold[peaks_by_scaffold$scaffold == scaf, ]
      pad <- pmax(0.1 * (scaf_peaks$Peak_End - scaf_peaks$Peak_Start), scaf_len * 0.015)
      scaf_peaks$band_start <- pmax(0, scaf_peaks$Peak_Start - pad)
      scaf_peaks$band_end   <- pmin(scaf_len, scaf_peaks$Peak_End + pad)

      jpeg(file.path(var_name,
                     paste0("Manhattan_40000_3_", scaf, "_", var_name, ".jpg")),
           width = 7, height = 3.5, units = "in", res = 300)
      print(
        ggplot(scaf_snps, aes(x = Position, y = -log10(p_final))) +
        geom_rect(data = scaf_peaks, inherit.aes = FALSE,
                  aes(xmin = band_start, xmax = band_end, ymin = -Inf, ymax = Inf),
                  fill = "grey80") +
        geom_point(size = 1.3, color = "black", alpha = 0.65) +
        geom_hline(yintercept = -log10(thresh_bonf), color = "darkblue",
                   linewidth = 1.1) +
        scale_x_continuous(breaks = x_breaks,
                           limits = c(0, scaf_len),
                           expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0.05)) +
        labs(x = NULL, y = NULL) +
        theme_bw(base_size = 11) +
        theme(legend.position  = "none",
              panel.grid.minor = element_blank(),
              axis.text.x      = element_blank(),
              axis.ticks.x     = element_blank(),
              axis.text.y      = element_blank(),
              axis.ticks.y     = element_blank())
      )
      dev.off()
      print(paste("Written per-scaffold Manhattan:", scaf))
    }
  }
}

}  # end Manhattan plots

# sink()
