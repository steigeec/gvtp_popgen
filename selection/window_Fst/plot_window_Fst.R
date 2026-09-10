#!/usr/bin/env Rscript

# Usage: Rscript plot_window_Fst.R [window_dir] [pair_table]
args <- commandArgs(trailingOnly = TRUE)
window_dir <- if (length(args) >= 1) args[1] else "/private/groups/meyerlab/esteiger/selection/window_Fst"
pair_table <- if (length(args) >= 2) args[2] else "/private/groups/meyerlab/esteiger/structure/Fst/populations/pairwise_auto/pairwise_Fst_populations.txt"

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")
library(data.table)
library(dplyr)
library(ggplot2)

FAI <- "/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked.fai"

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

# Genome-wide chromosome order/offsets and alternating grey/black coloring,
# matching the Manhattan plot convention in
# /private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma/find_outliers.R
# (scaffolds ordered by descending length, colors alternate per scaffold).
chroms <- fread(FAI, select = 1:2, col.names = c("chr", "len"))
chroms <- chroms %>% arrange(desc(len))
chroms$start <- cumsum(c(0, head(chroms$len, -1)))
REPS <- round(nrow(chroms) / 2) + 1
chrom_colors <- setNames(rep(c("grey36", "black"), REPS)[seq_len(nrow(chroms))], chroms$chr)

pairs <- read.table(pair_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

for (i in seq_len(nrow(pairs))) {
  p1 <- pairs$pop1[i]
  p2 <- pairs$pop2[i]
  infile <- file.path(window_dir, paste0(p1, ".", p2, "_window"))
  if (!file.exists(infile)) {
    message("Missing window file, skipping: ", infile)
    next
  }

  # realSFS fst stats2 writes a 4-name header (region chr midPos Nsites) for 5
  # actual data columns (region, chr, midPos, Nsites, Fst); fread would otherwise
  # shift every column name right by one. Skip the bad header and name explicitly.
  win <- fread(infile, skip = 1, col.names = c("region", "chr", "midPos", "Nsites", "Fst"))
  win <- win[is.finite(Fst) & chr %in% chroms$chr]
  if (nrow(win) == 0) {
    message("No usable windows, skipping: ", infile)
    next
  }

  win <- merge(win, chroms[, c("chr", "start")], by = "chr")
  win$pos <- win$midPos + win$start
  win$chr <- factor(win$chr, levels = chroms$chr)
  win <- win[order(win$pos)]

  outfile <- file.path(window_dir, paste0("window_Fst_", p1, ".", p2, ".jpg"))
  jpeg(outfile, width = 2000, height = 1000)
  print(ggplot(win, aes(x = pos, y = Fst, color = chr)) +
    geom_point(size = 2) +
    scale_color_manual(values = chrom_colors) +
    xlab("Position along the genome") + ylab("Fst") +
    ggtitle(paste(p1, "vs", p2)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0.03, 0.03)) +
    froggy_theme + theme(plot.title = element_text(size = 40, hjust = 0.5)))
  dev.off()
  message("Written: ", outfile)
}
