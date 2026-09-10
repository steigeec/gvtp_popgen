library(ggplot2)
library(dplyr)
library(ggrepel)

# Plots the reduced (rrvgo-collapsed) BP + MF GO enrichment results produced
# by GOenrich_goseq.R for a single gene-list scope. file_label must match
# GOenrich_goseq.R's own file_label for that run (gene_list_<file_label>.txt
# with the "gene_list_" prefix and ".txt" suffix stripped) — e.g.
# "fog_unfiltered", "BIO18_40000_3", "BIO15_40000_3_woInversion".
#
# Note: goseq's output has no "expected count" column, so there's no exact
# fold-enrichment here (that needs total DE-gene/universe counts that aren't
# in the per-term results table). Plotted instead: numDEInCat / numInCat,
# i.e. what fraction of a category's own genes are in the candidate set —
# still a fair measure of how concentrated a category's signal is.
#
# Usage: Rscript plot_GO_enrichment.R <file_label> [output_dir]
#   output_dir: directory holding GO_enrichment_results_<file_label>_*.tsv
#               (default: current directory)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript plot_GO_enrichment.R <file_label> [output_dir]")
}
file_label <- args[1]
output_dir <- if (length(args) >= 2) args[2] else "."

# ── Load data ──────────────────────────────────────────────────────────────────
bp_file <- file.path(output_dir, paste0("GO_enrichment_results_", file_label, "_BP_goseq_reduced.tsv"))
mf_file <- file.path(output_dir, paste0("GO_enrichment_results_", file_label, "_MF_goseq_reduced.tsv"))

bp <- if (file.exists(bp_file)) read.delim(bp_file, stringsAsFactors = FALSE) else NULL
mf <- if (file.exists(mf_file)) read.delim(mf_file, stringsAsFactors = FALSE) else NULL

if (is.null(bp) && is.null(mf)) {
    message("No BP or MF reduced enrichment results found for '", file_label,
            "' in ", output_dir, " — nothing to plot.")
    quit(save = "no", status = 0)
}

if (!is.null(bp)) bp$Ontology <- "BP"
if (!is.null(mf)) mf$Ontology <- "MF"
dat <- bind_rows(bp, mf)

dat <- dat %>%
  mutate(
    CategoryProportion = numDEInCat / numInCat,
    neglog10FDR         = -log10(FDR),
    Term                 = sub("\\.\\.\\.", "", term)   # strip trailing ellipsis if any
  )

dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)

# ── Plot 1: Dot plot ───────────────────────────────────────────────────────────
# Terms on y-axis ordered by category proportion within each ontology,
# x = category proportion, size = # candidate genes, color = -log10(FDR)

dat_ordered <- dat %>%
  arrange(Ontology, CategoryProportion) %>%
  mutate(Term = factor(Term, levels = unique(Term)))

p1 <- ggplot(dat_ordered, aes(x = CategoryProportion, y = Term,
                               size = numDEInCat, color = neglog10FDR)) +
  geom_point() +
  facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradient(low = "#b3cde3", high = "#1a1a8c",
                       name = expression(-log[10](FDR))) +
  scale_size_continuous(name = "Candidate\ngenes", range = c(4, 10)) +
  labs(
    title = paste0("GO Enrichment — ", file_label),
    x     = "Proportion of category genes in candidate set",
    y     = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text       = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey92"),
    legend.position  = "right"
  )

dotplot_pdf <- file.path(output_dir, "plots", paste0("GO_dotplot_", file_label, ".pdf"))
dotplot_png <- file.path(output_dir, "plots", paste0("GO_dotplot_", file_label, ".png"))
ggsave(dotplot_pdf, p1, width = 8, height = 5.5)
ggsave(dotplot_png, p1, width = 8, height = 5.5, dpi = 200)
message("Saved: ", dotplot_pdf, " / .png")


# ── Plot 2: Combined BP + MF bubble plot ──────────────────────────────────────
# x = category proportion, y = -log10(FDR), size = # candidate genes,
# color = ontology, label = term

p2 <- ggplot(dat, aes(x = CategoryProportion, y = neglog10FDR,
                       size = numDEInCat, color = Ontology, label = Term)) +
  geom_point(alpha = 0.75) +
  geom_text_repel(
    size          = 3.2,
    max.overlaps  = Inf,
    box.padding   = 0.4,
    show.legend   = FALSE
  ) +
  scale_color_manual(values = c(BP = "#2166ac", MF = "#d6604d"),
                     name = "Ontology") +
  scale_size_continuous(name = "Candidate\ngenes", range = c(4, 10)) +
  labs(
    title = paste0("GO Enrichment — ", file_label, " (BP + MF)"),
    x     = "Proportion of category genes in candidate set",
    y     = expression(-log[10](FDR))
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

bubble_pdf <- file.path(output_dir, "plots", paste0("GO_bubble_", file_label, ".pdf"))
bubble_png <- file.path(output_dir, "plots", paste0("GO_bubble_", file_label, ".png"))
ggsave(bubble_pdf, p2, width = 9, height = 6)
ggsave(bubble_png, p2, width = 9, height = 6, dpi = 200)
message("Saved: ", bubble_pdf, " / .png")
