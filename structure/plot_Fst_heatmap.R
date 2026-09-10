#!/usr/bin/env Rscript

## Pairwise population Fst heatmap (lower triangle only), sized for a small
## printed inset (~1x1 inch). Populations are marked by a colored dot along
## each axis (instead of text labels) to keep the grid legible at that scale.

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(ggplot2)

fst_file <- "/private/groups/meyerlab/esteiger/structure/Fst/populations/pairwise_auto/pairwise_Fst_populations.txt"
out_stem <- "/private/groups/meyerlab/esteiger/gvtp_popgen/structure/Fst_heatmap_populations"

fst_stat <- "Fst_weighted"   # switch to "Fst_unweighted" to plot the other estimator

# Pipeline population names (used in the Fst table) -> cluster names
name_map <- c(
  conception = "Point_Conception",
  cdfw       = "CDFW",
  arguello   = "VAFB_Arguello",
  strauss    = "Strauss",
  lionhead   = "VAFB_Lion_Head",
  Dangermond = "Dangermond"
)

cluster_levels <- c("VAFB_Lion_Head", "VAFB_Arguello", "Dangermond",
                     "Point_Conception", "Strauss", "CDFW")

# Bang Wong color palette
popcolors <- c("CDFW" = "#e69f00", "CDFW_MR" = "#e69f00", "Strauss"="#009e73",
               "Point_Conception"="#0072b2", "VAFB_Lion_Head"="#f0e442",
               "VAFB_Arguello"="#d55e00", "Dangermond"="#56b4e9", "Between"="black")

# ---------- Read and relabel ----------

fst_df <- read.delim(fst_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
fst_df$pop1 <- name_map[fst_df$pop1]
fst_df$pop2 <- name_map[fst_df$pop2]
if (any(is.na(fst_df$pop1)) || any(is.na(fst_df$pop2))) {
  stop("Unmapped population name(s) in Fst table -- check name_map.")
}

# ---------- Build the symmetric matrix, then keep only the lower triangle ----------

n <- length(cluster_levels)
mat <- matrix(NA_real_, n, n, dimnames = list(cluster_levels, cluster_levels))
for (i in seq_len(nrow(fst_df))) {
  p1 <- fst_df$pop1[i]; p2 <- fst_df$pop2[i]
  mat[p1, p2] <- fst_df[[fst_stat]][i]
  mat[p2, p1] <- fst_df[[fst_stat]][i]
}

# Numeric grid positions: col_i/row_i = 1..n in cluster_levels order.
# y is flipped (n+1-row_i) so the first population sits at the top.
plot_df <- expand.grid(row = cluster_levels, col = cluster_levels, stringsAsFactors = FALSE)
plot_df$col_i <- match(plot_df$col, cluster_levels)
plot_df$row_i <- match(plot_df$row, cluster_levels)
plot_df <- plot_df[plot_df$row_i > plot_df$col_i, ]   # lower triangle only, no self-comparisons
plot_df$fst <- mapply(function(r, c) mat[r, c], plot_df$row, plot_df$col)
plot_df$x <- plot_df$col_i
plot_df$y <- n + 1 - plot_df$row_i

# Axis dot strips: one dot per population that actually appears as a column /
# row in the lower triangle (drops CDFW from the x-axis and Lion's Head from
# the y-axis, since neither has any cells there).
cols_used <- sort(unique(plot_df$col_i))
rows_used <- sort(unique(plot_df$row_i))
col_dots <- data.frame(x = cols_used, y = 0, pop = cluster_levels[cols_used])
row_dots <- data.frame(x = 0, y = n + 1 - rows_used, pop = cluster_levels[rows_used])
col_dots$col <- popcolors[col_dots$pop]
row_dots$col <- popcolors[row_dots$pop]

fmt_fst <- function(x) sub("^0\\.", ".", sprintf("%.2f", x))   # ".05" not "0.05"

# White text on dark (high-Fst) tiles, black text on light tiles, so labels
# stay legible across the whole white-to-black gradient. Passed as a fixed
# vector (not via aes()) so it doesn't collide with the dots' colour scale.
mid <- mean(range(plot_df$fst))
plot_df$text_col <- ifelse(plot_df$fst > mid, "white", "black")

# Dot fill/outline colours are passed as fixed vectors (not through aes()) so
# they don't collide with the tiles' "fill" scale (continuous grey gradient).
# shape 21 gives a thin black outline around each dot -- needed since the
# palette's yellow (Lion's Head) is otherwise hard to see on white.
p <- ggplot() +
  geom_tile(data = plot_df, aes(x = x, y = y, fill = fst), color = "black", linewidth = 0.08) +
  geom_text(data = plot_df, aes(x = x, y = y, label = fmt_fst(fst)),
            colour = plot_df$text_col, size = 1.1) +
  geom_point(data = col_dots, aes(x = x, y = y),
             shape = 21, fill = col_dots$col, colour = "black", stroke = 0.2, size = 1.3) +
  geom_point(data = row_dots, aes(x = x, y = y),
             shape = 21, fill = row_dots$col, colour = "black", stroke = 0.2, size = 1.3) +
  scale_fill_gradient(low = "white", high = "black", guide = "none") +
  coord_fixed(xlim = c(-0.5, n + 0.5), ylim = c(-0.5, n + 0.5), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

# PDF (vector) is the primary output -- text stays crisp at 1x1"; the PNG is a
# quick-look raster at high dpi so small text doesn't blur.
ggsave(paste0(out_stem, ".pdf"), p, width = 1, height = 1, units = "in", bg = "white")
ggsave(paste0(out_stem, ".png"), p, width = 1, height = 1, units = "in", dpi = 600, bg = "white")

cat("Heatmap written to:\n  ", paste0(out_stem, ".pdf"), "\n  ", paste0(out_stem, ".png"), "\n")
