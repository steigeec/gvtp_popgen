#!/usr/bin/env Rscript

## Pairwise population Fst heatmap (lower triangle), sized for a ~1 x 1 inch
## inset. Populations are marked by coloured dots on the axes, not text labels.
## Dangermond cells are hatched because low coverage inflates its Fst
## (see mantel_IBD.R).

BASE_DIR <- "/private/groups/meyerlab/esteiger"
.libPaths(file.path(BASE_DIR, "software/R_libs/"))
library(ggplot2)

fst_file <- file.path(BASE_DIR, "structure/Fst/populations/pairwise_auto/pairwise_Fst_populations.txt")
out_stem <- file.path(BASE_DIR, "structure/Fst/populations/Fst_heatmap_populations")

fst_stat <- "Fst_weighted"   # switch to "Fst_unweighted" to plot the other estimator

# Fst table population names to cluster names
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

# Grid positions in cluster_levels order; y flipped so the first population is at the top.
plot_df <- expand.grid(row = cluster_levels, col = cluster_levels, stringsAsFactors = FALSE)
plot_df$col_i <- match(plot_df$col, cluster_levels)
plot_df$row_i <- match(plot_df$row, cluster_levels)
plot_df <- plot_df[plot_df$row_i > plot_df$col_i, ]   # lower triangle only, no self-comparisons
plot_df$fst <- mapply(function(r, c) mat[r, c], plot_df$row, plot_df$col)
plot_df$x <- plot_df$col_i
plot_df$y <- n + 1 - plot_df$row_i

# Axis dots only for populations with cells in the lower triangle
# (no CDFW on x, no Lion's Head on y).
cols_used <- sort(unique(plot_df$col_i))
rows_used <- sort(unique(plot_df$row_i))
col_dots <- data.frame(x = cols_used, y = 0, pop = cluster_levels[cols_used])
row_dots <- data.frame(x = 0, y = n + 1 - rows_used, pop = cluster_levels[rows_used])
col_dots$col <- popcolors[col_dots$pop]
row_dots$col <- popcolors[row_dots$pop]

# Hatch lines clipped to each tile, since ggplot2 has no pattern fill.
HATCH_POP <- "Dangermond"
HATCH_SPACING <- 0.25
hatch_cells <- plot_df[plot_df$row == HATCH_POP | plot_df$col == HATCH_POP, ]
hatch_df <- do.call(rbind, lapply(seq_len(nrow(hatch_cells)), function(i) {
  cx <- hatch_cells$x[i]; cy <- hatch_cells$y[i]
  # lines y - cy = (x - cx) - k across the tile [-0.5, 0.5]^2
  k <- seq(-1 + HATCH_SPACING, 1 - HATCH_SPACING, by = HATCH_SPACING)
  data.frame(x = cx + pmax(-0.5, -0.5 + k), xend = cx + pmin(0.5, 0.5 + k),
             y = cy + pmax(-0.5, -0.5 + k) - k, yend = cy + pmin(0.5, 0.5 + k) - k)
}))

fmt_fst <- function(x) sub("^0\\.", ".", sprintf("%.2f", x))   # ".05" not "0.05"

# White labels on dark tiles. Passed outside aes() so it doesn't clash with
# the dots' colour scale.
mid <- mean(range(plot_df$fst))
plot_df$text_col <- ifelse(plot_df$fst > mid, "white", "black")

# Dot colours also passed outside aes(), to keep them off the tile fill scale.
# Shape 21 adds a black outline so the yellow dot shows on white.
p <- ggplot() +
  geom_tile(data = plot_df, aes(x = x, y = y, fill = fst), color = "black", linewidth = 0.08) +
  geom_segment(data = hatch_df, aes(x = x, y = y, xend = xend, yend = yend),
               colour = "grey55", linewidth = 0.12) +
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

# PDF is the main output; PNG at 600 dpi for quick viewing.
ggsave(paste0(out_stem, ".pdf"), p, width = 1, height = 1, units = "in", bg = "white")
ggsave(paste0(out_stem, ".png"), p, width = 1, height = 1, units = "in", dpi = 600, bg = "white")

cat("Heatmap written to:\n  ", paste0(out_stem, ".pdf"), "\n  ", paste0(out_stem, ".png"), "\n")
