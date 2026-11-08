## Supplementary figure: salmon TPM of the GxE candidate genes across the 8
## tissue RNA-seq libraries. Plots the 51 genes in MANUSCRIPT_ORDER. Cell colour
## is each gene's log2(TPM + 1) scaled 0-1 across its 8 samples; the printed
## number is the TPM itself. Columns are in fixed tissue order (not clustered).
##
## Input:  candidate_gene_TPM_wide.tsv (aggregate_expression.py)
## Output: RNA_DIR/tissue_expression/results/supp_fig_tissue_heatmap.png
## R packages: ComplexHeatmap (Bioconductor), circlize (CRAN)
## Usage: Rscript supp_fig_tissue_heatmap.R

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

## ---- config ----------------------------------------------------------
RNA_DIR     <- "/scratch1/emma/RNA-seq"
WIDE_TSV    <- file.path(RNA_DIR, "tissue_expression/results/candidate_gene_TPM_wide.tsv")
GENE_BLOCKS <- 2     # number of side-by-side panels
OUT_PNG     <- file.path(RNA_DIR, "tissue_expression/results/supp_fig_tissue_heatmap.png")
FIG_WIDTH_IN  <- 6.5
FIG_HEIGHT_IN <- 4.3
PNG_RES     <- 400   # dpi

SAMPLES <- c("GVTP-1Root_S2","GVTP-2Root_S4","GVTP-5Root_S6","GVTPMRoot_S8",
             "GVTP-1Shoot_S1","GVTP-2Shoot_S3","GVTP-5Shoot_S5","GVTPMLeaf_S7")
TISSUE  <- c("Root","Root","Root","Root","Shoot","Shoot","Shoot","Leaf")
INDIV   <- c("1","2","5","M","1","2","5","M")

## ---- load + collapse to one row per gene ------------------------------
raw <- read.delim(WIDE_TSV, check.names = FALSE, stringsAsFactors = FALSE)
raw <- raw[!duplicated(raw$Gene_Name), ]

label <- ifelse(!is.na(raw$Acronym) & raw$Acronym != "" & raw$Acronym != "NA",
                raw$Acronym, raw$Gene_Name)

# Acronyms that are blank or out of date in gene_summary_all.tsv, checked against
# gene coordinates
acronym_overrides <- c(
  "IDmodified-gene-98221"  = "AP2/ERF-ANT",
  "g8215"                  = "ABCC",
  "g8291"                  = "RNASEH2B",
  "IDmodified-gene-287739" = "GGL2",
  "IDmodified-gene-21805"  = "CLPTM1",
  "IDmodified-gene-18840"  = "P450-CYP76",
  "IDmodified-gene-25607"  = "RhoGDI",
  "IDmodified-gene-23185"  = "CNIH4",
  "IDmodified-gene-115144" = "SEC1A"
)
hit <- match(raw$Gene_Name, names(acronym_overrides))
label[!is.na(hit)] <- acronym_overrides[hit[!is.na(hit)]]

# Number acronyms shared by more than one candidate (e.g. paralogous histones)
dup <- duplicated(label) | duplicated(label, fromLast = TRUE)
if (any(dup)) {
  ave_idx <- ave(seq_along(label), label, FUN = seq_along)
  label[dup] <- paste0(label[dup], " (", ave_idx[dup], ")")
}
stopifnot(!any(duplicated(label)))

tpm <- as.matrix(raw[, SAMPLES])
rownames(tpm) <- label

## TGD2 (fog candidate) is missing from gene_summary_all.tsv; its TPM is taken
## from salmon's quant.sf for the matching gene, DEINVI-g8012
## (DEINVI_00017:56439469-56448202).
tgd2_tpm <- c(6.020517, 8.469273, 6.803258, 7.008907, 9.219180, 13.623607, 11.405188, 10.330727)
names(tgd2_tpm) <- SAMPLES
tpm <- rbind(tpm, "TGD2" = tgd2_tpm[SAMPLES])

log2tpm <- log2(tpm + 1)

## ---- row-normalize for color; keep raw TPM for the printed values -----
row_sd <- apply(log2tpm, 1, sd)
row_min <- apply(log2tpm, 1, min)
row_max <- apply(log2tpm, 1, max)

row_scaled <- (log2tpm - row_min) / (row_max - row_min)
row_scaled[row_sd == 0, ] <- 0  # flat or undetected genes shown as "low"

## ---- genes plotted, in figure order -------------------------------------
MANUSCRIPT_ORDER <- c(
  "PAT11","NPF8.4","FPGS3","GLUD1","SCPL13","ANL2","CSA1","CYP21","ADF5",
  "L16","L14","S11","rpoA","RTFL6","NST1","SYP32","UCH2","AP2/ERF-ANT",
  "HBT","TGD2","ABCC","RNASEH2B","ANN2","APK3","GFS12","KIN12C","APC6",
  "ALG2","PAP18","GLO1","SCPL17","GGL2","LIP1","GNS1","SNF1","CLPTM1",
  "PSI2","MAM33","P450-CYP76","RhoGDI","DRP3B","EID1","AGL24","SAC51",
  "RPOT3","SAPK2","AHL10","CNIH4","SEC1A","PAT14","NRS/ER"
)
stopifnot(all(MANUSCRIPT_ORDER %in% rownames(log2tpm)))
ordered_genes <- MANUSCRIPT_ORDER

## ---- split into side-by-side blocks ------------------------------------
n <- length(ordered_genes)
block_size <- ceiling(n / GENE_BLOCKS)
block_id <- ceiling(seq_len(n) / block_size)
gene_blocks <- split(ordered_genes, block_id)

## ---- per-row color scale (0-1, min-max within each gene) ---------------
col_fun <- colorRamp2(
  c(0, 0.25, 0.5, 0.75, 1),
  c("#eaf1fb", "#9ec5f4", "#5598e7", "#2a78d6", "#0d366b")
)

tissue_cols <- c(Root = "#935c30", Shoot = "#3c7a49", Leaf = "#7a9a34")

## ---- build one Heatmap per block ----------------------------------------
## Blocks hold different genes, so they are placed with grid viewports rather
## than ComplexHeatmap's `ht1 + ht2` (which needs shared rows).
fmt_tpm <- function(x) ifelse(x < 10, sprintf("%.1f", x), sprintf("%.0f", round(x)))

make_block <- function(genes) {
  disp <- row_scaled[genes, , drop = FALSE]  # drives cell color
  raw  <- tpm[genes, , drop = FALSE]         # drives printed text
  Heatmap(
    disp,
    name              = "rel_expr",
    col               = col_fun,
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    show_row_dend     = FALSE,
    row_labels        = rownames(disp),
    row_names_gp      = gpar(fontsize = 6),
    row_names_side    = "left",
    column_split      = factor(TISSUE, levels = c("Root", "Shoot", "Leaf")),
    column_title_gp   = gpar(fontsize = 7, fontface = "bold"),
    column_gap        = unit(0.8, "mm"),
    column_labels     = INDIV,
    column_names_gp   = gpar(fontsize = 6),
    column_names_rot  = 0,
    column_names_centered = TRUE,
    top_annotation    = HeatmapAnnotation(
      tissue = TISSUE,
      col = list(tissue = tissue_cols),
      show_annotation_name = FALSE,
      show_legend = FALSE,
      simple_anno_size = unit(1.6, "mm")
    ),
    rect_gp           = gpar(col = "white", lwd = 0.4),
    border            = TRUE,
    show_heatmap_legend = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      ink <- if (disp[i, j] > 0.55) "white" else "#14171a"
      grid.text(fmt_tpm(raw[i, j]), x, y, gp = gpar(fontsize = 4.6, col = ink))
    }
  )
}

blocks <- lapply(gene_blocks, make_block)

## Colour is relative within each gene, so the legend shows Low/High only.
legend <- Legend(
  col_fun = col_fun,
  title = "Expression\n(row-scaled)",
  at = c(0, 1),
  labels = c("Low", "High"),
  title_gp = gpar(fontsize = 7, fontface = "plain"),
  labels_gp = gpar(fontsize = 6.5),
  legend_height = unit(14, "mm")
)

## ---- render ----------------------------------------------------------
png(OUT_PNG, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN, units = "in", res = PNG_RES)
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow = 1, ncol = GENE_BLOCKS + 1,
  widths = unit(c(rep(1, GENE_BLOCKS), 0.4), "null")
)))
for (i in seq_along(blocks)) {
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
  draw(blocks[[i]], newpage = FALSE)
  popViewport()
}
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = GENE_BLOCKS + 1))
draw(legend, just = "left")
popViewport()
popViewport()
dev.off()
message("Wrote ", OUT_PNG)
