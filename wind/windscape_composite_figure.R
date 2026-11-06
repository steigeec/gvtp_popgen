## 2x2 composite figure for the windscape connectivity analysis:
## A = population-scale net-flow map (HRRR), B = population Fst vs wind connectivity,
## C = Tranquillon Mtn net-flow map (100 m WindNinja), D = individual Fst vs wind connectivity.
##
## Inputs (from windscape_connectivity.R): hrrr_wind_rose.tif,
##   windninja_tranquillon_wind_rose.tif, pop_hours_hrrr.txt,
##   strauss_hours_windninja.txt; OSM coastline (build_coastline.sh);
##   population and individual pairwise Fst; field_data.txt
## Outputs: windscape/composite_panels/panel_[a-d].jpg,
##   windscape/windscape_connectivity_composite.jpg
## R libraries: R_libs_windscape, R_libs (windscape, terra, jpeg)
## Usage: Rscript windscape_composite_figure.R
##
## Scatter plots are unadjusted, so they include isolation by distance.
## Panels are rendered separately and composited with rasterImage(), because
## terra's plot() does not work with par(mfrow).

BASE_DIR <- "/private/groups/meyerlab/esteiger"
.libPaths(c(file.path(BASE_DIR, "software/R_libs_windscape"),
            file.path(BASE_DIR, "software/R_libs/")))
library(windscape)
library(terra)
library(jpeg)

WIND_DIR <- file.path(BASE_DIR, "wind")
WS_DIR <- file.path(WIND_DIR, "windscape")
PANEL_DIR <- file.path(WS_DIR, "composite_panels")
dir.create(PANEL_DIR, showWarnings = FALSE)

# ---------- Shared styling ----------
PANEL_W <- 6; PANEL_H <- 5.4; PANEL_RES <- 300; PANEL_PS <- 12
PANEL_W_NARROW <- PANEL_W * 0.6  # B and D are 60% of the width of A and C
CEX_AXIS <- 1.4
CEX_MTEXT <- 1.3
CEX_LABEL <- 0.9
FLOW_COL <- hcl.colors(140, "YlGnBu", rev = TRUE)[1:100]
XLAB_CONN <- "Wind connectivity\nln(mean flow, 1/h)"

name_map <- c(
  "Point_Conception" = "conception", "CDFW" = "cdfw",
  "VAFB_Arguello"    = "arguello",   "Strauss" = "strauss",
  "VAFB_Lion_Head"   = "lionhead",   "Dangermond" = "Dangermond"
)
pop_colors <- c(
  cdfw = "#e69f00", strauss = "#009e73", conception = "#0072b2",
  lionhead = "#f0e442", arguello = "#d55e00", Dangermond = "#56b4e9"
)
pop_labels <- c(
  lionhead = "Lion's Head", arguello = "Pt. Arguello", strauss = "Tranquillon",
  Dangermond = "Dangermond", conception = "Pt. Conception", cdfw = "Gaviota"
)

# The 13 excluded sample IDs have untrustworthy coordinates.
BAD_COORD_IDS <- c("4548", "4549", "4550", "4551", "4552", "4553", "4554",
                    "4555", "4556", "4557", "4559", "4561", "4562")
field <- read.delim(file.path(BASE_DIR, "field_data.txt"), header = TRUE, sep = "\t")
field$SampleID <- as.character(field$SampleID)
field_clean <- field[!field$SampleID %in% BAD_COORD_IDS, ]
centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field_clean, FUN = mean)
centroids$pop <- name_map[centroids$Clusters]
coords <- as.matrix(centroids[, c("Long", "Lat")])
rownames(coords) <- centroids$pop

read_hours <- function(file) {
  h <- as.matrix(read.delim(file.path(WS_DIR, file), row.names = 1, check.names = FALSE))
  off <- row(h) != col(h)
  h[off & h == 0] <- min(h[off & h > 0])  # coincident individuals
  h
}
connectivity <- function(hours) log(pairwise_means(1 / hours))

# Net-flow map. Arrows at every `thin`-th cell; `arrow_scale` is degrees per m/s.
# Overlays are drawn inside fun= because terra resets par when plot() returns.
flow_map <- function(rose_file, display_extent, arrow_scale, thin, arrow_lwd, path, add_points,
                     left_mar = 5) {
  nf <- net_flow(wind_rose(file.path(WS_DIR, rose_file)))
  rate <- sqrt(nf[["u"]]^2 + nf[["v"]]^2)
  arr <- as.data.frame(c(nf, rate), xy = TRUE, na.rm = TRUE)
  names(arr)[5] <- "rate"
  arr <- arr[arr$x >= display_extent[1] & arr$x <= display_extent[2] &
             arr$y >= display_extent[3] & arr$y <= display_extent[4], ]
  col_i <- round((arr$x - min(arr$x)) / res(nf)[1])
  row_i <- round((arr$y - min(arr$y)) / res(nf)[2])
  arr <- arr[col_i %% thin == 0 & row_i %% thin == 0, ]

  jpeg(path, width = PANEL_W, height = PANEL_H, units = "in", res = PANEL_RES, pointsize = PANEL_PS)
  plot(crop(rate, display_extent, snap = "out"), ext = display_extent, col = FLOW_COL, asp = NA,
       mar = c(2.2, left_mar, 2.5, 5.5),
       pax = list(cex.axis = CEX_AXIS, las = 1, side = c(1, 2)),
       plg = list(title = "Net wind flow (m/s)", title.cex = CEX_MTEXT * 0.85, cex = CEX_AXIS * 0.75),
       fun = function() {
         arrows(arr$x, arr$y, arr$x + arr$u * arrow_scale, arr$y + arr$v * arrow_scale,
                length = 0.05, lwd = arrow_lwd, col = "black")
         add_points()
       })
  dev.off()
}

# `pad` keeps large symbols from being clipped; pad = 0 fits the x axis to the data.
scatter_frame <- function(x, y, pad = 0.08) {
  par(mar = c(6.5, 4.6, 1, 1))
  plot(x, y, type = "n", xlim = extendrange(x, f = pad), ylim = extendrange(y, f = 0.08),
       xaxs = if (pad == 0) "i" else "r",
       xlab = "", ylab = "", cex.axis = CEX_AXIS, cex.lab = CEX_MTEXT)
  title(ylab = expression(F[ST]/(1-F[ST])), line = 2.8, cex.lab = CEX_MTEXT)
  title(xlab = XLAB_CONN, line = 4.8, cex.lab = CEX_MTEXT)
}

# ---------- Panel A: population-scale windscape map (HRRR) ----------

label_pos <- setNames(rep(3, nrow(coords)), rownames(coords))
label_pos["arguello"] <- 1

coastline <- vect(file.path(BASE_DIR, "env/coastline/osm_coastline_study_area.gpkg"))

panel_a_path <- file.path(PANEL_DIR, "panel_a.jpg")
flow_map("hrrr_wind_rose.tif", ext(-120.7, -120.1, 34.4, 34.9), arrow_scale = 0.006, thin = 2,
         arrow_lwd = 1.2, path = panel_a_path, add_points = function() {
  lines(coastline, col = "black", lwd = 1.5)
  points(coords[, "Long"], coords[, "Lat"], pch = 21, bg = pop_colors[rownames(coords)], col = "black", cex = 1.6, lwd = 1)
  text(coords[, "Long"], coords[, "Lat"], labels = pop_labels[rownames(coords)],
       pos = label_pos[rownames(coords)], cex = CEX_LABEL, font = 2)
})

# ---------- Panel B: population Fst vs windscape wind connectivity (HRRR) ----------

fst_df <- read.delim(
  file.path(BASE_DIR, "structure/Fst/populations/pairwise_auto/pairwise_Fst_populations.txt"),
  header = TRUE, sep = "\t")
pops <- sort(unique(c(fst_df$pop1, fst_df$pop2)))
fst_mat <- matrix(0, length(pops), length(pops), dimnames = list(pops, pops))
fst_mat[cbind(fst_df$pop1, fst_df$pop2)] <- fst_df$Fst_weighted
fst_mat[cbind(fst_df$pop2, fst_df$pop1)] <- fst_df$Fst_weighted
fst_lin <- fst_mat / (1 - fst_mat)
conn <- connectivity(read_hours("pop_hours_hrrr.txt")[pops, pops])

pairs_idx <- which(upper.tri(fst_lin), arr.ind = TRUE)
pair_df <- data.frame(pop1 = pops[pairs_idx[, 1]], pop2 = pops[pairs_idx[, 2]],
                      conn = conn[pairs_idx], fst = fst_lin[pairs_idx])

panel_b_path <- file.path(PANEL_DIR, "panel_b.jpg")
jpeg(panel_b_path, width = PANEL_W_NARROW, height = PANEL_H, units = "in", res = PANEL_RES, pointsize = PANEL_PS)
scatter_frame(pair_df$conn, pair_df$fst)
# Dangermond pairs are hollow (Fst inflated by low coverage; excluded from the tests).
dang <- pair_df$pop1 == "Dangermond" | pair_df$pop2 == "Dangermond"
partner <- ifelse(pair_df$pop1 == "Dangermond", pair_df$pop2, pair_df$pop1)
points(pair_df$conn[!dang], pair_df$fst[!dang], pch = 21, cex = 2.6, lwd = 4.5,
       bg = pop_colors[pair_df$pop1[!dang]], col = pop_colors[pair_df$pop2[!dang]])
points(pair_df$conn[dang], pair_df$fst[dang], pch = 21, cex = 2.6, lwd = 4.5,
       bg = "white", col = pop_colors[partner[dang]])
dev.off()

# ---------- Panel C: Tranquillon Mtn windscape map (WindNinja) ----------

strauss <- field[field$Clusters == "Strauss", ]
pad_deg <- 0.003
panel_c_path <- file.path(PANEL_DIR, "panel_c.jpg")
flow_map("windninja_tranquillon_wind_rose.tif",
         ext(min(strauss$Long) - pad_deg, max(strauss$Long) + pad_deg,
             min(strauss$Lat) - pad_deg, max(strauss$Lat) + pad_deg),
         arrow_scale = 0.0004, thin = 3, arrow_lwd = 1.5, path = panel_c_path, left_mar = 6.5, add_points = function() {
  points(strauss$Long, strauss$Lat, pch = 21, bg = "white", col = "black", cex = 0.7, lwd = 0.6)
})

# ---------- Panel D: individual Fst vs windscape wind connectivity (WindNinja) ----------

bam_ids <- sub("_rmdup\\.bam$", "", readLines(file.path(BASE_DIR, "mapping/merged/bamlist_strauss.txt")))
fst_all <- read.delim(file.path(BASE_DIR, "structure/Fst/individuals/pairwise_Fst.txt"), header = TRUE)
fst_all$sample1 <- as.character(fst_all$sample1); fst_all$sample2 <- as.character(fst_all$sample2)
fst_s <- fst_all[fst_all$sample1 %in% bam_ids & fst_all$sample2 %in% bam_ids, ]
fst_s <- fst_s[fst_s$sample1 != "6432" & fst_s$sample2 != "6432", ]  # duplicate-sampled plant
ids <- sort(unique(c(fst_s$sample1, fst_s$sample2)))
fst_mat_ind <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
fst_mat_ind[cbind(fst_s$sample1, fst_s$sample2)] <- fst_s$Fst_weighted
fst_mat_ind[cbind(fst_s$sample2, fst_s$sample1)] <- fst_s$Fst_weighted
fst_lin_ind <- fst_mat_ind / (1 - fst_mat_ind)
conn_ind <- connectivity(read_hours("strauss_hours_windninja.txt")[ids, ids])

x_ind <- conn_ind[lower.tri(conn_ind)]
y_ind <- fst_lin_ind[lower.tri(fst_lin_ind)]

panel_d_path <- file.path(PANEL_DIR, "panel_d.jpg")
jpeg(panel_d_path, width = PANEL_W_NARROW, height = PANEL_H, units = "in", res = PANEL_RES, pointsize = PANEL_PS)
scatter_frame(x_ind, y_ind, pad = 0)
points(x_ind, y_ind, pch = 19, col = adjustcolor("#009e73", alpha.f = 0.35), cex = 0.6)
abline(lm(y_ind ~ x_ind), col = "black", lwd = 2)
dev.off()

# ---------- Composite: 2x2 grid with panel labels ----------

TOTAL_W <- PANEL_W + PANEL_W_NARROW
TOTAL_H <- PANEL_H * 2
col_split <- PANEL_W / TOTAL_W

out_path <- file.path(WS_DIR, "windscape_connectivity_composite.jpg")
jpeg(out_path, width = TOTAL_W, height = TOTAL_H, units = "in", res = PANEL_RES)
panel_info <- list(
  list(path = panel_a_path, letter = "A", fig = c(0, col_split, 0.5, 1)),
  list(path = panel_b_path, letter = "B", fig = c(col_split, 1, 0.5, 1)),
  list(path = panel_c_path, letter = "C", fig = c(0, col_split, 0, 0.5)),
  list(path = panel_d_path, letter = "D", fig = c(col_split, 1, 0, 0.5))
)
for (i in seq_along(panel_info)) {
  pinfo <- panel_info[[i]]
  par(fig = pinfo$fig, mar = c(0, 0, 0, 0), new = (i > 1))
  img <- readJPEG(pinfo$path)
  plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  rasterImage(img, 0, 0, 1, 1)
  text(0.015, 0.99, pinfo$letter, cex = 2.2, font = 2, adj = c(0, 1), xpd = NA)
}
dev.off()
cat("Saved composite ->", out_path, "\n")
