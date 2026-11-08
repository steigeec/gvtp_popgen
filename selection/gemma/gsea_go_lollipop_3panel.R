#!/usr/bin/env Rscript --vanilla

# Lollipop plot of nominally enriched GO terms (gsea_go_stricter_nominal.R output)
# for fog (A), precipitation seasonality / BIO15 (B) and summer precipitation /
# BIO18 (C), side by side on one page. Length = -log10 p; colour = functional
# theme (PANEL_THEMES, by GO ID; unthemed terms grey).
#
# Usage: Rscript gsea_go_lollipop_3panel.R [tag] [top_n]
#   tag: "" (default) or e.g. "sunflower", matching gsea_go_stricter_nominal.R
# Writes: GO_lollipop_3panel<_tag>.pdf / .png in the gemma directory

BASE_DIR <- "/private/groups/meyerlab/esteiger"

suppressMessages({
    .libPaths(file.path(BASE_DIR, "software/R_libs"))
    library(dplyr); library(ggplot2); library(patchwork); library(stringr)
})

GEMMA_DIR  <- file.path(BASE_DIR, "selection/gemma")
PANELS     <- c(fog = "A  Fog", BIO15 = "B  Precipitation seasonality", BIO18 = "C  Summer precipitation")
# Functional themes per panel (by GO ID); interpretive groupings of nominal
# terms, assigned only where the term name points clearly to the theme
PANEL_THEMES <- list(
    fog = list(
        "Cell wall and carbohydrate remodeling"  = c("GO:0016759", "GO:0000272", "GO:0031501", "GO:0015925"),
        "Altered light perception"               = c("GO:0090227", "GO:0009584"),
        "Osmoregulation and leaf surface traits" = c("GO:0010325", "GO:1904278"),
        "Immune function"                        = c("GO:0051502"),
        "Developmental timing"                   = c("GO:0010331", "GO:0045962"),
        "Lipid and energy metabolism"            = c("GO:0010896", "GO:0004806", "GO:0072350")),
    BIO15 = list(
        "Cell wall and carbohydrate remodeling"  = c("GO:0009251", "GO:0016759", "GO:0016160", "GO:0008878",
                                                     "GO:1902503", "GO:0006112"),
        "Chloroplast function and photoprotection" = c("GO:0045037", "GO:0031897", "GO:0080183", "GO:0008465"),
        "Immune function"                        = c("GO:0051502", "GO:0036202"),
        "Redox and detoxification"               = c("GO:0010731"),
        "Developmental timing"                   = c("GO:1900056")),
    BIO18 = list(
        "Chloroplast function and photoprotection" = c("GO:0009654", "GO:0009543", "GO:0033614", "GO:0009767"),
        "Immune function"                        = c("GO:0050687", "GO:0031332", "GO:0046741", "GO:0045548"),
        "Redox and detoxification"               = c("GO:1990748", "GO:0046185", "GO:1901654"),
        "Developmental timing"                   = c("GO:0045544", "GO:0048575"),
        "Osmoregulation and leaf surface traits" = c("GO:1904276")))
# one colour per theme, shared across panels (Okabe-Ito; yellow darkened for text)
THEME_COLORS <- c("Cell wall and carbohydrate remodeling"    = "#E69F00",
                  "Altered light perception"                 = "#56B4E9",
                  "Osmoregulation and leaf surface traits"   = "#009E73",
                  "Immune function"                          = "#D55E00",
                  "Developmental timing"                     = "#CC79A7",
                  "Lipid and energy metabolism"              = "#0072B2",
                  "Chloroplast function and photoprotection" = "#9E8B00",
                  "Redox and detoxification"                 = "#000000",
                  "Other"                                    = "grey70")
OTHER_TEXT <- "grey50"     # label colour for unthemed terms (black is a theme colour)
PAGE_W_IN  <- 8.5      # at 183 mm print width, text scales to ~85%
PAGE_H_IN  <- 7.9
TEXT_PT    <- 6.5      # term-name font size (journal minimum is 5-7 pt)
WRAP_CHARS <- 48
ABBREVIATE <- TRUE
SHOW_GO_ID <- TRUE     # "name (GO:0000000)"; the ID is kept on one line

args  <- commandArgs(trailingOnly = TRUE)
tag   <- if (length(args) >= 1) args[1] else ""
TOP_N <- if (length(args) >= 2) as.numeric(args[2]) else Inf
sfx   <- if (tag != "") paste0("_", tag) else ""

abbrev <- function(x) {
    if (!ABBREVIATE) return(x)
    x <- str_replace_all(x, "positive regulation of", "pos. reg. of")
    x <- str_replace_all(x, "negative regulation of", "neg. reg. of")
    x <- str_replace_all(x, "regulation of", "reg. of")
    x <- str_replace_all(x, "biosynthetic process", "biosynthesis")
    x <- str_replace_all(x, "catabolic process", "catabolism")
    x <- str_replace_all(x, "metabolic process", "metabolism")
    x
}

dat <- lapply(names(PANELS), function(v) {
    d <- read.table(file.path(GEMMA_DIR, v, paste0("GSEA_GO_stricter_nominal_", v, sfx, ".tsv")),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
    d <- d[order(d$pval), ]
    if (is.finite(TOP_N)) d <- head(d, TOP_N)
    d$var <- v
    d
})

panel <- function(d, title, themes) {
    lab <- abbrev(d$term)
    # wrap the name, then append the ID as one unbreakable token
    if (SHOW_GO_ID) lab <- paste0(lab, "\u00a0(", sub("GO:", "GO:", d$pathway), ")")
    theme_of <- setNames(rep(names(themes), lengths(themes)), unlist(themes))
    d$theme <- factor(ifelse(d$pathway %in% names(theme_of), theme_of[d$pathway], "Other"),
                      levels = names(THEME_COLORS))
    d <- d %>% mutate(neglog10p = -log10(pval), label = str_wrap(lab, WRAP_CHARS)) %>%
        arrange(desc(neglog10p)) %>% mutate(label = factor(label, levels = rev(unique(label))))
    # label colours in axis (factor-level) order
    lab_col <- with(d[match(levels(d$label), d$label), ],
                    ifelse(theme == "Other", OTHER_TEXT, THEME_COLORS[as.character(theme)]))
    ggplot(d, aes(x = neglog10p, y = label)) +
        geom_segment(aes(x = 0, xend = neglog10p, yend = label, color = theme), linewidth = 0.55) +
        geom_point(aes(color = theme), size = 1.4) +
        # legend lists this panel's themes only (in THEME_COLORS order), not "Other"
        scale_color_manual(values = THEME_COLORS, name = NULL,
                           breaks = intersect(names(THEME_COLORS), names(themes))) +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 1.8))) +
        scale_x_continuous(expand = expansion(mult = c(0.01, 0.06))) +
        labs(x = expression(-log[10](p)), y = NULL, title = title) +
        theme_bw(base_size = 7) +
        theme(axis.text.y = element_text(size = TEXT_PT, lineheight = 0.8, color = lab_col),
              axis.text.x = element_text(size = 6),
              plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 3)),
              plot.title.position = "plot",
              panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
              legend.position = "bottom", legend.justification = "right",
              legend.text = element_text(size = 6.5), legend.key.height = unit(8, "pt"),
              legend.key.spacing.y = unit(0, "pt"), legend.margin = margin(2, 0, 0, 0),
              plot.margin = margin(2, 4, 2, 2))
}

p <- wrap_plots(mapply(panel, dat, PANELS, PANEL_THEMES[names(PANELS)], SIMPLIFY = FALSE), nrow = 1)

out <- file.path(GEMMA_DIR, paste0("GO_lollipop_3panel", sfx))
ggsave(paste0(out, ".pdf"), p, width = PAGE_W_IN, height = PAGE_H_IN)
ggsave(paste0(out, ".png"), p, width = PAGE_W_IN, height = PAGE_H_IN, dpi = 300)
cat("Saved:", paste0(out, ".pdf / .png"), sprintf("(%s terms; %.1f x %.1f in)\n",
    paste(sapply(dat, nrow), collapse = " / "), PAGE_W_IN, PAGE_H_IN))
