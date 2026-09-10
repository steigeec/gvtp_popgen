#!/usr/bin/env Rscript --vanilla

# Lollipop plot of gsea_go_stricter_nominal.R's stricter-threshold, rrvgo-
# reduced GO terms (replaces an earlier gene-overlap network figure for the
# same table -- with this few terms and this little leading-edge gene
# overlap between them at p < 0.005, the network mostly showed isolated
# points rather than real structure, so a ranked lollipop reads better).
#
# Deliberately minimal: no title/subtitle, no gene-name annotations.
# Grayscale by ontology (BP black, MF dark grey, CC light grey) rather than
# the earlier hue-based palette. A bottom-anchored legend covers both
# encodings (ontology color, term-size point scale) without competing with
# the term-name labels for horizontal space.
#
# No FDR-significant terms exist in this dataset -- like
# gsea_go_stricter_nominal.R, this is nominal/hypothesis-generating, not a
# substitute for the FDR-based result.
#
# Usage: Rscript gsea_go_lollipop.R <var_name> [output_dir] [top_n] [tag] [width_in] [height_in]
#   e.g. Rscript gsea_go_lollipop.R fog fog 12 "" 2 3
# Expects <output_dir>/GSEA_GO_stricter_nominal_<var_name>.tsv
# (gsea_go_stricter_nominal.R output).
# Writes: <output_dir>/plots/GO_lollipop_<var_name>.pdf / .png

suppressMessages({
    library(dplyr)
    library(ggplot2)
})

ONT_COLORS <- c(BP = "black", MF = "grey40", CC = "grey75")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript gsea_go_lollipop.R <var_name> [output_dir] [top_n] [tag] [width_in] [height_in]")
}
var_name   <- args[1]
output_dir <- if (length(args) >= 2) args[2] else var_name
TOP_N      <- if (length(args) >= 3) as.numeric(args[3]) else Inf
# tag must match gsea_go_stricter_nominal.R's tag for this run (e.g.
# "sunflower"), so this reads/writes the matching tagged file set.
tag        <- if (length(args) >= 4) args[4] else ""
tag_suffix <- if (tag != "") paste0("_", tag) else ""
WIDTH_IN   <- if (length(args) >= 5) as.numeric(args[5]) else NA  # NA -> auto (see plot_height/plot_width below)
HEIGHT_IN  <- if (length(args) >= 6) as.numeric(args[6]) else NA

infile <- file.path(output_dir, paste0("GSEA_GO_stricter_nominal_", var_name, tag_suffix, ".tsv"))
if (!file.exists(infile)) {
    stop("Missing ", infile, " -- run gsea_go_stricter_nominal.R for '", var_name, "' first.")
}
dat <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
dat <- dat[order(dat$pval), ]
if (is.finite(TOP_N) && nrow(dat) > TOP_N) {
    print(paste0("Keeping top ", TOP_N, " of ", nrow(dat), " terms by p-value"))
    dat <- dat[seq_len(TOP_N), ]
}

dat$neglog10p <- -log10(dat$pval)
dat <- dat %>% arrange(desc(neglog10p)) %>%
    mutate(term = factor(term, levels = rev(unique(term))))

p <- ggplot(dat, aes(x = neglog10p, y = term)) +
    geom_segment(aes(x = 0, xend = neglog10p, y = term, yend = term, color = ontology),
                linewidth = 0.7) +
    geom_point(aes(color = ontology, size = size)) +
    scale_color_manual(values = ONT_COLORS, name = "Ontology") +
    scale_size_continuous(name = "Term size (genes)", range = c(1.2, 4)) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
    labs(x = expression(-log[10](p)), y = NULL) +
    theme_bw(base_size = 15) +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin = margin(4, 14, 4, 4),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(t = 2)
    )

dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
lolli_pdf <- file.path(output_dir, "plots", paste0("GO_lollipop_", var_name, tag_suffix, ".pdf"))
lolli_png <- file.path(output_dir, "plots", paste0("GO_lollipop_", var_name, tag_suffix, ".png"))
plot_width  <- if (!is.na(WIDTH_IN))  WIDTH_IN  else 11
plot_height <- if (!is.na(HEIGHT_IN)) HEIGHT_IN else max(2, nrow(dat) * 0.22 + 0.5) + 0.9  # +0.9in for the bottom legend
ggsave(lolli_pdf, p, width = plot_width, height = plot_height, limitsize = FALSE)
ggsave(lolli_png, p, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
print(paste0("Saved: ", lolli_pdf, " / .png (", nrow(dat), " terms, ",
             plot_width, "in x ", plot_height, "in)"))
