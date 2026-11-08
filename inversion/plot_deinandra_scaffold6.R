#!/usr/bin/env Rscript
# Synteny plots for Deinandra scaffold 6 (DEINVI_00006) from the finished
# GENESPACE run:
#  1. a riparian plot with scaffold 6 as the anchor and the 5 other genomes as tracks
#  2. pairwise dot plots of scaffold 6 vs each other genome, from the
#     syntenicHits/*.synHits.txt.gz files (GENESPACE's own dotplots cannot be
#     subset to one scaffold)
# Output: genespace/custom_plots/*.jpg
# Usage: Rscript plot_deinandra_scaffold6.R
# Conda env: genespace
#   built by annotation/setup_genespace.sh: conda env with r-base=4.3,
#   MCScanX and GENESPACE's Bioconductor/CRAN dependencies, then
#   devtools::install_github("jtlovell/GENESPACE") (GENESPACE 1.3.1)

BASE_DIR <- "/private/groups/meyerlab/esteiger"

library(GENESPACE)
library(data.table)
library(ggplot2)

wd <- file.path(BASE_DIR, "inversions/synteny/genespace")
outDir <- file.path(wd, "custom_plots")
scaffold <- "DEINVI_00006"

# Putative inversion in bp; plots use gene rank order on the x axis, so
# membership is computed from the bp columns and used only for colouring.
invStart <- 43e6
invEnd <- 50e6

# display names only; genomeIDs are join keys and stay unchanged
prettyName <- function(x) gsub("_", " ", x)

genomeIDs <- c(
  "Deinandra", "Helianthus_annuus", "Helianthus_anomalus",
  "Lactuca_sativa", "Arabidopsis_thaliana", "Solanum_lycopersicum")

gpar <- init_genespace(
  wd = wd,
  genomeIDs = genomeIDs,
  ploidy = rep(1, length(genomeIDs)),
  path2mcscanx = "/private/home/esteiger/miniconda3/envs/genespace/bin",
  rawOrthofinderDir = file.path(wd, "orthofinder"),
  nCores = 1)

##############################################################################
# 1. Riparian plot, scaffold 6 only: inversion in GENESPACE's default
#    highlight colour, rest of the scaffold grey. Earlier highlightBed rows are
#    drawn on top, so the inversion row must come first.
cat("Building combined riparian plot for", scaffold, "...\n")

hb <- data.frame(
  genome = c("Deinandra", "Deinandra"),
  chr = c(scaffold, scaffold),
  start = c(invStart, 0),
  end = c(invEnd, Inf),
  color = c("#C4645C", "grey70"))

rip <- plot_riparian(
  gsParam = gpar,
  genomeIDs = genomeIDs,
  refGenome = "Deinandra",
  highlightBed = hb,
  backgroundColor = NULL,
  pdfFile = NULL)

ripPlot <- rip$plotData$ggplotObj

# genome names are the y-scale tick labels, so relabel that scale in place
yScaleIdx <- which(vapply(ripPlot$scales$scales, function(s) "y" %in% s$aesthetics, logical(1)))
ripPlot$scales$scales[[yScaleIdx]]$labels <- prettyName(ripPlot$scales$scales[[yScaleIdx]]$labels)

ggsave(
  filename = file.path(outDir, "Deinandra_scaffold6_riparian.jpg"),
  plot = ripPlot,
  device = "jpeg", width = 8, height = length(genomeIDs) * 1.4,
  dpi = 300, limitsize = FALSE)
cat("  wrote Deinandra_scaffold6_riparian.jpg\n")

##############################################################################
# 2. Pairwise dot plots, scaffold 6 vs each other genome
synDir <- gpar$paths$syntenicHits
otherGenomes <- setdiff(genomeIDs, "Deinandra")

for(og in otherGenomes){
  cat("Building dot plot: Deinandra scaffold 6 vs", og, "...\n")

  f1 <- file.path(synDir, sprintf("Deinandra_vs_%s.synHits.txt.gz", og))
  f2 <- file.path(synDir, sprintf("%s_vs_Deinandra.synHits.txt.gz", og))
  f <- if(file.exists(f1)) f1 else f2
  if(!file.exists(f))
    stop("could not find synHits file for Deinandra vs ", og, " (checked ", f1, " and ", f2, ")")

  hits <- read_synHits(f)

  # Deinandra can be genome1 or genome2; standardize so it is always the x axis
  if(hits$genome1[1] == "Deinandra"){
    hits <- subset(hits, chr1 == scaffold)
    setnames(hits,
      c("chr1","start1","end1","ord1","chr2","start2","end2","ord2"),
      c("chrD","startD","endD","ordD","chrO","startO","endO","ordO"))
  } else {
    hits <- subset(hits, chr2 == scaffold)
    setnames(hits,
      c("chr2","start2","end2","ord2","chr1","start1","end1","ord1"),
      c("chrD","startD","endD","ordD","chrO","startO","endO","ordO"))
  }

  if(nrow(hits) == 0){
    cat("  ** no syntenic hits between Deinandra scaffold 6 and", og, "-- skipping plot\n")
    next
  }

  # drop target chromosomes with < 10 hits (scattered noise; their facets
  # would be illegibly thin), as GENESPACE's minGenes2plot does
  nPerChr <- hits[,.N, by = chrO]
  keepChrs <- nPerChr$chrO[nPerChr$N >= 10]
  nDropped <- uniqueN(hits$chrO) - length(keepChrs)
  hits <- subset(hits, chrO %in% keepChrs)
  if(nDropped > 0)
    cat("  (dropped", nDropped, "target chr(s) with <10 hits to avoid squished facets)\n")

  # colour by overlap of the Deinandra gene with the inversion window
  hits[,inInversion := pmax(startD, invStart) <= pmin(endD, invEnd)]

  p <- ggplot(hits, aes(x = ordD, y = ordO, color = inInversion)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "grey70"),
                        guide = "none") +
    facet_grid(~ chrO, scales = "free_x", space = "free_x") +
    labs(
      x = sprintf("Deinandra %s: gene rank order", scaffold),
      y = sprintf("%s: gene rank order", prettyName(og)),
      title = sprintf("Deinandra %s vs %s", scaffold, prettyName(og))) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(size = 7), panel.grid.minor = element_blank())

  ggsave(
    filename = file.path(outDir, sprintf("Deinandra_scaffold6_vs_%s.jpg", og)),
    plot = p, device = "jpeg", width = 9, height = 6, dpi = 300)
  cat("  wrote Deinandra_scaffold6_vs_", og, ".jpg\n", sep = "")
}

cat("Done. All plots under:", outDir, "\n")
