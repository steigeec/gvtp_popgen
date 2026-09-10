library(ggplot2)
library(dplyr)

# Renders a styled summary table of the reduced (rrvgo-collapsed) BP + MF GO
# enrichment results produced by GOenrich_goseq.R, for a single gene-list
# scope. file_label must match GOenrich_goseq.R's own file_label for that
# run (gene_list_<file_label>.txt with the "gene_list_" prefix and ".txt"
# suffix stripped) — e.g. "fog_unfiltered", "BIO18_40000_3".
#
# Note: unlike a much earlier hardcoded draft of this script, all data below
# is read from the actual enrichment output — there is no baked-in gene list
# or GO term data here.
#
# Usage: Rscript plot_GO_table.R <file_label> [output_dir]
#   output_dir: directory holding GO_enrichment_results_<file_label>_*.tsv
#               (default: current directory)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript plot_GO_table.R <file_label> [output_dir]")
}
file_label <- args[1]
output_dir <- if (length(args) >= 2) args[2] else "."

dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)

# ── Load data ──────────────────────────────────────────────────────────────────
bp_file <- file.path(output_dir, paste0("GO_enrichment_results_", file_label, "_BP_goseq_reduced.tsv"))
mf_file <- file.path(output_dir, paste0("GO_enrichment_results_", file_label, "_MF_goseq_reduced.tsv"))

bp <- if (file.exists(bp_file)) read.delim(bp_file, stringsAsFactors = FALSE) else NULL
mf <- if (file.exists(mf_file)) read.delim(mf_file, stringsAsFactors = FALSE) else NULL

if (is.null(bp) && is.null(mf)) {
    message("No BP or MF reduced enrichment results found for '", file_label,
            "' in ", output_dir, " — nothing to tabulate.")
    quit(save = "no", status = 0)
}

if (!is.null(bp)) bp$Ontology <- "BP"
if (!is.null(mf)) mf$Ontology <- "MF"
dat <- bind_rows(bp, mf) %>% arrange(FDR)

dat <- dat %>%
  mutate(
    Term    = sub("\\.\\.\\.", "", term),
    FoldEnr = numDEInCat / numInCat
  )

dat$FDR_fmt <- formatC(dat$FDR,    format = "e", digits = 2)
dat$FE_fmt  <- sprintf("%.2f",     dat$FoldEnr)

# ── Layout: columns and rows ───────────────────────────────────────────────────
cols <- data.frame(
  col  = c("Ont","GO_ID","Term","FE","FDR","Genes"),
  xmin = c(0,    0.9,    2.4,   9.4, 10.5, 12.0),
  xmax = c(0.9,  2.4,    9.4,   10.5,12.0, 20.0),
  stringsAsFactors = FALSE
)
cols$xmid <- (cols$xmin + cols$xmax) / 2

n <- nrow(dat)
header_y <- n + 1

rows_df <- lapply(seq_len(n), function(i) {
  row <- dat[i, ]
  is_even <- i %% 2 == 0
  row_bg <- ifelse(is_even, "#f4f4f4", "white")
  ont_bg <- ifelse(row$Ontology == "BP", "#d6e8f7", "#fde4db")

  data.frame(
    row_i  = i,
    y      = n + 1 - i,
    col    = cols$col,
    xmin   = cols$xmin,
    xmax   = cols$xmax,
    xmid   = cols$xmid,
    fill   = c(ont_bg, rep(row_bg, 5)),
    label  = c(row$Ontology, row$category, row$Term,
               row$FE_fmt, row$FDR_fmt, row$Genes),
    hjust  = c(0.5, 0, 0, 0.5, 0.5, 0),
    xpos   = c(cols$xmid[1],
               cols$xmin[2] + 0.15,
               cols$xmin[3] + 0.15,
               cols$xmid[4],
               cols$xmid[5],
               cols$xmin[6] + 0.15),
    fontface = c("bold", rep("plain", 5)),
    color    = c(ifelse(row$Ontology == "BP", "#1a5a9a", "#c0392b"),
                 rep("grey15", 5)),
    stringsAsFactors = FALSE
  )
})
rows_df <- do.call(rbind, rows_df)

hdr_labels <- c("Ont.", "GO ID", "GO Term", "Prop.\nin cat.", "FDR", "Genes")
hdr_df <- data.frame(
  y     = header_y,
  xmin  = cols$xmin,
  xmax  = cols$xmax,
  xmid  = cols$xmid,
  fill  = "#2c4770",
  label = hdr_labels,
  hjust = c(0.5, 0, 0, 0.5, 0.5, 0),
  xpos  = c(cols$xmid[1],
             cols$xmin[2] + 0.15,
             cols$xmin[3] + 0.15,
             cols$xmid[4],
             cols$xmid[5],
             cols$xmin[6] + 0.15),
  color = "white",
  fontface = "bold",
  stringsAsFactors = FALSE
)

# ── Plot ───────────────────────────────────────────────────────────────────────
p <- ggplot() +
  geom_rect(data = rows_df,
            aes(xmin = xmin, xmax = xmax, ymin = y - 0.5, ymax = y + 0.5, fill = fill),
            color = "grey80", linewidth = 0.25) +
  scale_fill_identity() +
  geom_rect(data = hdr_df,
            aes(xmin = xmin, xmax = xmax,
                ymin = y - 0.55, ymax = y + 0.55, fill = fill),
            color = NA) +
  geom_text(data = rows_df,
            aes(x = xpos, y = y, label = label,
                hjust = hjust, color = color, fontface = fontface),
            size = 2.6, lineheight = 0.9) +
  scale_color_identity() +
  geom_text(data = hdr_df,
            aes(x = xpos, y = y, label = label, hjust = hjust),
            color = "white", fontface = "bold", size = 2.9, lineheight = 0.85) +
  annotate("text", x = 0, y = 0.1,
           label = "Prop. in cat. = candidate genes / total genes in category.  FDR = Benjamini-Hochberg corrected p-value.",
           hjust = 0, size = 2.2, color = "grey55") +
  scale_x_continuous(limits = c(0, 20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.1, header_y + 0.6), expand = c(0, 0)) +
  labs(title = paste0("GO Enrichment Summary — ", file_label)) +
  theme_void(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", size = 13, margin = margin(b = 6)),
    plot.margin  = margin(10, 10, 6, 10)
  )

table_pdf <- file.path(output_dir, "plots", paste0("GO_table_", file_label, ".pdf"))
table_png <- file.path(output_dir, "plots", paste0("GO_table_", file_label, ".png"))
ggsave(table_pdf, p, width = 14, height = max(3, n * 0.4 + 1.5))
ggsave(table_png, p, width = 14, height = max(3, n * 0.4 + 1.5), dpi = 220)
message("Saved: ", table_pdf, " / .png")
