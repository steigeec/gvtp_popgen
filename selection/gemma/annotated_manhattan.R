#!/usr/bin/env Rscript
#
# Manhattan plots of the 40000_3 retained peaks (same Bonferroni and window filtering as
# find_outliers.R), labelling every candidate gene in <var_name>/gene_summary_<var_name>.tsv
# that SnpEff assigned to a retained SNP, with its acronym or else ortholog name (genes with
# neither are unlabelled). Labels sit near the top of the plot with a segment to the peak.
#
# The inversion region (DEINVI_00006:43-50Mb) is excluded from the whole-
# genome plot's labels and shown only on the DEINVI_00006 scaffold-level
# plot.
#
# Usage: Rscript annotated_manhattan.R <pheno_num> <var_name>  (run from the gemma directory)

BASE_DIR <- "/private/groups/meyerlab/esteiger"

args <- commandArgs(trailingOnly=TRUE)
pheno_num <- args[1]
var_name  <- args[2]

dir.create(var_name, showWarnings=FALSE, recursive=TRUE)

.libPaths(file.path(BASE_DIR, "software/R_libs"))
library(bit64)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

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

dic <- read.table("scaffold_snp_map.txt")
dic <- unique(dic[,c("V1", "V3")])
scaf2chr <- setNames(as.character(dic$V1), dic$V3)
chr2scaf <- setNames(dic$V3, as.character(dic$V1))

gemma <- read.table(paste0("output/pheno", pheno_num, "_lmm_results.assoc.txt"), header=T)
names(gemma)[1:3] <- c("Chromosome", "SNP", "Position")
print(paste("No rows in input:", nrow(gemma)))

# Genomic control, as in find_outliers.R
chisq_obs <- qchisq(gemma$p_lrt, df=1, lower.tail=FALSE)
lambda <- median(chisq_obs, na.rm=TRUE) / qchisq(0.5, df=1)
print(paste("Lambda:", round(lambda, 3)))

if (lambda > 1.05) {
  gemma$p_final <- pchisq(chisq_obs / lambda, df=1, lower.tail=FALSE)
  print("Applied genomic control")
} else {
  gemma$p_final <- gemma$p_lrt
  print("No genomic control needed")
}

thresh_bonf <- 0.05 / nrow(gemma)
bonf_hits <- sum(gemma$p_final < thresh_bonf, na.rm=TRUE)
print(paste("Bonferroni hits:", bonf_hits))

fdr_results <- p.adjust(gemma$p_final, method = "fdr")
fdr_hits <- sum(fdr_results <= 0.05, na.rm=TRUE)
print(paste("FDR hits:", fdr_hits))

if (bonf_hits == 0 && fdr_hits == 0) {
  print("No outliers from Bonferroni or FDR — skipping Manhattan plot")
  quit(save="no", status=0)
}

gemma$out <- "no"
gemma$out[gemma$p_final < thresh_bonf] <- "out"

data <- gemma
data_no <- data[data$out=="no",]
data_no <- data_no[sample(nrow(data_no), nrow(data_no)/2), ]

pks <- as.data.table(data[data$out=="out",])

if (nrow(pks) == 0) {
  print("No Bonferroni outliers — nothing to plot")
  quit(save="no", status=0)
}

# Sliding-window retention, as in find_outliers.R
win <- 40000; step <- 20000
win_table <- list()
for (ch in unique(pks$Chromosome)) {
  chrlen <- max(pks[Chromosome==ch, Position])
  starts <- seq(1, chrlen, by = step)
  ends <- pmin(starts + win - 1, chrlen)
  dtw <- data.table(Chromosome = ch, start = starts, end = ends)
  setkey(dtw, start, end)
  win_table[[as.character(ch)]] <- dtw
}
wins <- rbindlist(win_table)
setkey(pks, Position, Position)
wins2 <- wins[, .(Chromosome, wstart = start, wend = end)]
setnames(wins2, c("Chromosome","wstart","wend"))
setkey(wins2, wstart, wend)
win_counts <- list()
for (ch in unique(wins2$Chromosome)) {
  w <- wins2[Chromosome == ch]
  s <- pks[Chromosome == ch]
  if (nrow(s) == 0) next
  setnames(s, "Position", "pos")
  s[, pos2 := pos]
  setkey(s, pos, pos2)
  setkey(w, wstart, wend)
  ol <- foverlaps(w, s, by.x = c("wstart", "wend"), by.y = c("pos", "pos2"), type = "any", nomatch = 0L)
  cnt <- ol[, .N, by = .(Chromosome, wstart, wend)]
  setnames(cnt, c("Chromosome", "start", "end", "N"))
  win_counts[[ch]] <- cnt
  setnames(s, "pos", "Position")
  s[, pos2 := NULL]
}
win_counts_dt <- rbindlist(win_counts)
win_kept <- win_counts_dt[N >= 3]
print(paste("Windows with >=2 SNPs within 40kb:", nrow(win_kept)))

setnames(win_kept, c("start", "end"), c("wstart", "wend"), skip_absent = TRUE)

pks[, pos2 := Position]
setkey(pks, Chromosome, Position, pos2)
setkey(win_kept, Chromosome, wstart, wend)

retained <- unique(foverlaps(pks, win_kept, by.x = c("Chromosome", "Position", "pos2"),
                      by.y = c("Chromosome", "wstart", "wend"), nomatch = 0L))
retained[, c("wstart", "wend", "N", "pos2") := NULL]
retained <- unique(retained)
pks[, pos2 := NULL]
print(paste("Retained SNPs (multi-SNP peaks):", nrow(retained)))

if (nrow(retained) == 0) {
  print("No filtered (40000_3) peaks — nothing to plot")
  quit(save="no", status=0)
}
retained$out <- "out"

# Add cumulative genomic position columns, as in find_outliers.R
add_positions <- function(man) {
  chroms <- man %>%
    group_by(Chromosome) %>%
    summarise(max_pos = max(Position)) %>%
    arrange(desc(max_pos))
  chroms$start <- 0
  lng <- 0
  for (ROW in seq_len(nrow(chroms))) {
    chroms$start[ROW] <- lng
    lng <- lng + chroms$max_pos[ROW]
  }
  man$Chromosome   <- as.character(man$Chromosome)
  chroms$Chromosome <- as.character(chroms$Chromosome)
  man <- merge(man, chroms[, c("Chromosome", "start")], by = "Chromosome")
  man$pos        <- man$Position + man$start
  man$Chromosome <- factor(man$Chromosome, levels = chroms$Chromosome)
  man <- man %>% arrange(pos)
  list(man = man, chroms = chroms)
}

res_ret <- add_positions(rbind(data_no, as.data.frame(retained)))
man    <- res_ret$man
chroms <- res_ret$chroms

# ── Gene name/acronym lookup ───────────────────────────────────────────────
# Label = Acronym, else cleaned ortholog name; genes with neither are dropped.
gene_summary_file <- file.path(var_name, paste0("gene_summary_", var_name, ".tsv"))
if (file.exists(gene_summary_file)) {
  gene_summary <- read.table(gene_summary_file, header=TRUE, sep="\t",
                              quote="", comment.char="", stringsAsFactors=FALSE)
  # All candidate genes; inclusion is decided by top_snp_per_gene (retained SNPs only),
  # not the In_40000_3_Peak flag.
  genes <- gene_summary
  genes$Gene_Coordinates <- genes$Coordinates
  # Bare AGI locus tags (e.g. "At5g09550") are not used as labels
  is_agi_locus_id <- function(x) grepl("^At[0-9CM]G[0-9]+$", x, ignore.case = TRUE)

  # Reject Acronym values without any alphanumeric character (e.g. "; ; ;")
  has_acronym  <- !is.na(genes$Acronym)   & genes$Acronym   != "" & genes$Acronym   != "NA" &
                   grepl("[A-Za-z0-9]", genes$Acronym) & !is_agi_locus_id(genes$Acronym)

  # Orthologs can be "; "-joined; generic names ("AT2G25670 protein", "Uncharacterized
  # protein ...") are dropped part by part. GENERIC_ORTHOLOG_EXACT: manually reviewed
  # eggNOG strings that name no gene.
  GENERIC_ORTHOLOG_EXACT <- c(
    "floral organ abscission",
    "Mitochondrial protein",
    "Belongs to the TRAFAC class myosin-kinesin ATPase superfamily. Kinesin family",
    "Belongs to the adaptor complexes small subunit family",
    "Belongs to the peptidase C1 family"
  )
  is_generic_ortholog <- function(x) {
    grepl("^AT[0-9CM]G[0-9]+(\\s+protein)?$", x, ignore.case = TRUE) |
    grepl("uncharacterized protein", x, ignore.case = TRUE) |
    is_agi_locus_id(x) |
    x %in% GENERIC_ORTHOLOG_EXACT
  }
  # Use only the first ortholog name, truncated, so labels stay inside the panel
  MAX_LABEL_CHARS <- 45
  clean_orthologs <- function(x) {
    if (is.na(x) || x == "" || x == "NA") return(NA_character_)
    parts <- trimws(strsplit(x, ";")[[1]])
    parts <- parts[nchar(parts) > 0 & !is_generic_ortholog(parts)]
    if (length(parts) == 0) return(NA_character_)
    best <- parts[1]
    if (nchar(best) > MAX_LABEL_CHARS) best <- paste0(substr(best, 1, MAX_LABEL_CHARS - 1), "…")
    best
  }
  genes$Orthologs_clean <- vapply(genes$Orthologs, clean_orthologs, character(1))
  has_ortholog <- !is.na(genes$Orthologs_clean)

  genes$Label <- ifelse(has_acronym, genes$Acronym, ifelse(has_ortholog, genes$Orthologs_clean, NA))
  genes <- genes[!is.na(genes$Label), ]
} else {
  print(paste(gene_summary_file, "not found — peaks will be left unlabeled"))
  genes <- data.frame(Gene_ID=character(), Gene_Coordinates=character(), Gene_Name=character(), Label=character())
}
print(paste("Candidate genes for", var_name, "with a usable label (pre retained-SNP matching):", nrow(genes)))

# Assign labels to height tiers alternating by x order (via nudge_y), so x-adjacent labels
# start at different heights; ggrepel then adjusts locally.
assign_tier_nudge <- function(lab_df, y_floor, y_ceiling, y_col, n_tiers = 3) {
  n <- nrow(lab_df)
  if (n == 0) return(numeric(0))
  n_tiers <- min(n_tiers, n)
  tier <- (rank(lab_df$pos, ties.method = "first") - 1) %% n_tiers
  tier_height <- if (n_tiers > 1) (y_ceiling - y_floor) / (n_tiers - 1) else 0
  target_y <- y_floor + tier * tier_height
  target_y - lab_df[[y_col]]
}

# For each gene, the retained SNP with the lowest p that SnpEff assigned to it
# (snp_to_gene_<var_name>.tsv from build_snp_gene_map.py; most outlier SNPs are
# intergenic, so coordinate overlap would miss most genes).
snp_gene_file <- file.path(var_name, paste0("snp_to_gene_", var_name, ".tsv"))
snp_gene <- read.table(snp_gene_file, header=TRUE, sep="\t", quote="", comment.char="",
                        stringsAsFactors=FALSE)
snp_gene_map <- new.env()
for (i in seq_len(nrow(snp_gene))) {
  key <- paste(snp_gene$Scaffold[i], snp_gene$Position[i])
  assign(key, strsplit(snp_gene$Gene_IDs[i], ";")[[1]], envir = snp_gene_map)
}
genes_for_snp <- function(scaffold, position) {
  get0(paste(scaffold, position), envir = snp_gene_map, ifnotfound = character(0))
}

top_snp_per_gene <- function(genes, snps) {
  # Match on Gene_ID (Gene_Name is a display label and can differ)
  target_ids <- genes$Gene_ID
  # gene_id -> best (lowest-p) retained-SNP row seen so far
  best <- new.env()
  for (i in seq_len(nrow(snps))) {
    scaf <- chr2scaf[[as.character(snps$Chromosome[i])]]
    if (is.null(scaf)) next
    gids <- genes_for_snp(scaf, snps$Position[i])
    gids <- gids[gids %in% target_ids]
    if (length(gids) == 0) next
    for (gid in gids) {
      cur <- get0(gid, envir = best, ifnotfound = NULL)
      if (is.null(cur) || snps$p_final[i] < cur$p_final) {
        assign(gid, snps[i, c("Chromosome", "Position", "p_final")], envir = best)
      }
    }
  }
  rows <- list()
  for (i in seq_len(nrow(genes))) {
    gid <- genes$Gene_ID[i]
    top <- get0(gid, envir = best, ifnotfound = NULL)
    if (is.null(top)) next
    rows[[length(rows) + 1]] <- data.frame(
      Chromosome  = top$Chromosome,
      Position    = top$Position,
      p_final     = top$p_final,
      Gene_Name   = genes$Gene_Name[i],
      Label       = genes$Label[i],
      has_acronym = !is.na(genes$Acronym[i]) & genes$Acronym[i] != "" & genes$Acronym[i] != "NA" &
                    grepl("[A-Za-z0-9]", genes$Acronym[i]) & !is_agi_locus_id(genes$Acronym[i])
    )
  }
  if (length(rows) > 0) do.call(rbind, rows) else
    data.frame(Chromosome=character(), Position=numeric(), p_final=numeric(),
               Gene_Name=character(), Label=character(), has_acronym=logical())
}

retained_df <- as.data.frame(retained)
candidates  <- top_snp_per_gene(genes, retained_df)

# Putative inversion region (excluded from whole-genome plot labels)
INV_SCAFFOLD <- "DEINVI_00006"
INV_START <- 43000000
INV_END   <- 50000000
inv_chr <- scaf2chr[[INV_SCAFFOLD]]

# ── Labels: one per candidate gene (a peak can carry several) ──────────────
label_cols <- c("Chromosome", "Position", "p_final", "Gene_Name", "Label", "has_acronym")
labels <- if (nrow(candidates) > 0) unique(candidates[, label_cols]) else candidates[, label_cols]

# Report of labelled genes that use an ortholog name instead of an acronym
needs_acronym <- labels[!labels$has_acronym, c("Gene_Name", "Label")]
write.table(needs_acronym[order(needs_acronym$Gene_Name), ],
            file.path(var_name, paste0("labeled_genes_needing_acronym_", var_name, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(paste("Labeled genes without a real acronym:", nrow(needs_acronym)))

if (nrow(labels) > 0) {
  labels$Chromosome <- as.character(labels$Chromosome)
  labels <- merge(labels, chroms[, c("Chromosome", "start")], by = "Chromosome")
  labels$pos <- labels$Position + labels$start
}
print(paste("Genes labeled:", nrow(labels)))

# Labels inside the inversion window are shown only on the DEINVI_00006 scaffold plot
main_labels <- labels
if (!is.null(inv_chr) && nrow(labels) > 0) {
  in_inversion <- as.character(labels$Chromosome) == inv_chr &
                  labels$Position >= INV_START & labels$Position <= INV_END
  main_labels <- labels[!in_inversion, ]
}

# ── Plot: filtered (40000_3) peaks, gene labels pinned near the top ───────
REPS <- round(length(levels(man$Chromosome)) / 2) + 1
y_top <- max(-log10(man$p_final), na.rm = TRUE)
# geom_text_repel's ylim does not extend the axis, so set the panel limit explicitly,
# above the label ceiling to leave room for the text itself.
y_label_ceiling <- y_top * 2.0
y_panel_ceiling <- y_top * 4.0

jpeg(file.path(var_name, paste0("Manhattan-geneLabel_40000_3_", var_name, ".jpg")), width = 2000, height = 1000)
p <- ggplot(man) +
  labs(x = NULL) + ylab("-log10(p-value)") +
  geom_point(aes(x = pos, y = -log10(p_final), col = Chromosome, size = out, shape = out)) +
  scale_color_manual(values = rep(c("grey36", "black"), REPS)) +
  geom_hline(yintercept = -log10(thresh_bonf), col = "grey50", lwd = 3) +
  scale_size_manual(values  = c("out" = 5, "no" = 2)) +
  scale_shape_manual(values = c("out" = 19, "no" = 19)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.03, 0.03), limits = c(0, y_panel_ceiling)) +
  froggy_theme

if (nrow(main_labels) > 0) {
  # Labels near the top, staggered over 8 tiers; direction = "y" keeps each at its peak's x
  main_labels$logp <- -log10(main_labels$p_final)
  main_nudge <- assign_tier_nudge(main_labels, y_top * 1.05, y_label_ceiling, "logp", n_tiers = 8)
  p <- p + geom_text_repel(data = main_labels,
                            aes(x = pos, y = logp, label = Label),
                            size = 4, fontface = "italic", color = "black",
                            direction = "y",
                            nudge_y = main_nudge,
                            ylim = c(y_top * 1.05, y_label_ceiling),
                            box.padding = 0.3, point.padding = 0.2, force = 3,
                            segment.color = "grey30", segment.size = 0.5,
                            max.overlaps = Inf, min.segment.length = 0)
}
print(p)
dev.off()

# ── Inversion-scaffold Manhattan, gene labels pinned near the top ─────────
# Whole DEINVI_00006 scaffold in the per-scaffold style of find_outliers.R, labelled from
# the full `labels` table (including inversion-window genes), with larger y-axis text,
# a heavier threshold line and a flatter aspect ratio.
if (!is.null(inv_chr)) {
  inv_retained <- retained_df[as.character(retained_df$Chromosome) == inv_chr, ]

  if (nrow(inv_retained) > 0) {
    inv_snps   <- data[as.character(data$Chromosome) == inv_chr, ]
    inv_labels <- labels[as.character(labels$Chromosome) == inv_chr, ]

    FAI <- file.path(BASE_DIR, "genome/DEINVI_genome.fasta.masked.fai")
    fai_df <- read.table(FAI, header = FALSE,
                         col.names = c("scaffold", "length", "offset", "linebases", "linewidth"),
                         stringsAsFactors = FALSE)
    scaf_len <- fai_df$length[fai_df$scaffold == INV_SCAFFOLD]
    x_breaks <- seq(0, scaf_len, by = 10e6)

    # Label and panel ceilings, as in the whole-genome plot
    inv_y_top           <- max(-log10(inv_snps$p_final), na.rm = TRUE)
    inv_y_label_ceiling <- inv_y_top * 2.0
    inv_y_panel_ceiling <- inv_y_top * 4.0

    jpeg(file.path(var_name, paste0("Manhattan-geneLabel_40000_3_", INV_SCAFFOLD, "_", var_name, ".jpg")),
         width = 8, height = 2.3, units = "in", res = 300)
    p_inv <- ggplot(inv_snps, aes(x = Position, y = -log10(p_final))) +
      geom_point(size = 0.5, color = "black", alpha = 0.65) +
      geom_hline(yintercept = -log10(thresh_bonf), color = "grey50",
                 linetype = "dashed", linewidth = 1.5) +
      scale_x_continuous(breaks = x_breaks,
                         limits = c(0, scaf_len),
                         expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0.05), limits = c(0, inv_y_panel_ceiling)) +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 11) +
      theme(legend.position  = "none",
            panel.grid.minor = element_blank(),
            axis.text.x      = element_blank(),
            axis.ticks.x     = element_blank(),
            axis.text.y      = element_text(size = 17.6))

    if (nrow(inv_labels) > 0) {
      inv_labels$logp <- -log10(inv_labels$p_final)
      inv_nudge <- assign_tier_nudge(inv_labels, inv_y_top * 1.05, inv_y_label_ceiling, "logp", n_tiers = 8)
      p_inv <- p_inv + geom_text_repel(data = inv_labels,
                                        aes(x = Position, y = logp, label = Label),
                                        size = 2, fontface = "italic", color = "black",
                                        direction = "y",
                                        nudge_y = inv_nudge,
                                        ylim = c(inv_y_top * 1.05, inv_y_label_ceiling),
                                        box.padding = 0.3, point.padding = 0.15, force = 3,
                                        segment.color = "grey30", segment.size = 0.3,
                                        max.overlaps = Inf, min.segment.length = 0)
    }
    print(p_inv)
    dev.off()
    print(paste("Written inversion-scaffold Manhattan, labeled genes:", nrow(inv_labels)))
  } else {
    print("No retained (40000_3) peaks on the inversion scaffold — skipping inversion Manhattan plot")
  }
}

# ── Plot-ready per-scaffold data for the stacked inversion figures ─────────
# (gvtp_popgen/inversion/DEINVI_0006_stacked.py): same SNPs, threshold and labels as the
# per-scaffold plot above (non-outliers plus retained 40000_3 SNPs, not thinned).
EXPORT_SCAFFOLDS <- c("DEINVI_00003", "DEINVI_00006", "DEINVI_00008", "DEINVI_00017")
fai_all <- read.table(file.path(BASE_DIR, "genome/DEINVI_genome.fasta.masked.fai"),
                      header = FALSE, stringsAsFactors = FALSE)[, 1:2]
retained_df <- as.data.frame(retained)
for (scaf in EXPORT_SCAFFOLDS) {
  if (!scaf %in% names(scaf2chr)) next
  chr  <- scaf2chr[[scaf]]
  snps <- data[as.character(data$Chromosome) == chr, ]
  ret  <- retained_df[as.character(retained_df$Chromosome) == chr, ]
  keep <- snps$out == "no" | snps$Position %in% ret$Position
  fwrite(data.table(Position = snps$Position[keep], p_final = snps$p_final[keep],
                    thresh_bonf = thresh_bonf, scaf_len = fai_all$V2[fai_all$V1 == scaf]),
         file.path(var_name, paste0("manhattan_", scaf, "_", var_name, ".tsv.gz")), sep = "\t")
  scaf_labels <- if (nrow(labels) > 0) labels[as.character(labels$Chromosome) == chr, ] else labels
  fwrite(as.data.table(scaf_labels[, c("Position", "p_final", "Label")]),
         file.path(var_name, paste0("manhattan_labels_", scaf, "_", var_name, ".tsv")), sep = "\t")
  print(paste("Exported", scaf, ":", sum(keep), "SNPs,", nrow(scaf_labels), "labels"))
}

print(paste("Done:", var_name))
