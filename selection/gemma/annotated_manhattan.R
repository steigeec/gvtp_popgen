#!/usr/bin/env Rscript
#
# Re-runs the same Bonferroni / 40000_3 peak-filtering logic as
# find_outliers.R, then labels EVERY gene flagged In_40000_3_Peak == "Yes"
# in <var_name>/gene_summary_<var_name>.tsv (i.e. every LD-curated
# candidate gene actually tied to a retained peak for this variable, not
# a subset chosen by adaptive-role/GO-term evidence or by picking one
# winner per peak) with its short acronym, falling back to its ortholog
# name when no acronym is resolved; genes with neither are left
# unlabeled. All labels are plain black -- no GO-term color-coding or
# legend (dropped per explicit request: the plots are dense enough with
# every retained-peak gene shown that a color scheme added more noise
# than signal, and the legend also happened to trigger a genuine ggrepel/
# grid rendering bug the one time it had real terms to draw — see git
# history for the crash this replaced). Label text is deliberately small
# since a peak can now carry several labels rather than just one. Labels
# are pinned near the top of the plot (small vertical stagger to avoid
# overlap) with a thin segment pointing down to the actual peak, rather
# than sitting on top of the peak itself.
#
# The inversion region (DEINVI_00006:43-50Mb) is excluded from the whole-
# genome plot's labels and shown only on the DEINVI_00006 scaffold-level
# plot.

args <- commandArgs(trailingOnly=TRUE)
pheno_num <- args[1]
var_name  <- args[2]

dir.create(var_name, showWarnings=FALSE, recursive=TRUE)

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")
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

# Genomic control — identical to find_outliers.R
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

# Sliding-window retention: keep SNPs with >= 2 neighbours within 40kb — identical to find_outliers.R
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

# Add cumulative genomic position columns — identical to find_outliers.R
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
# Per-variable sheet (built by compile_gene_summary.py from the UNFILTERED
# LD-curated candidate gene list), filtered to In_40000_3_Peak == "Yes" --
# every gene actually tied to a retained (Bonferroni + 40kb-window-
# density-filtered) peak for this variable, not a subset chosen by
# adaptive-role/GO-term evidence. Label = Acronym when resolved, else the
# Ortholog name when present; genes with neither are excluded entirely
# (no raw-gene-id fallback).
gene_summary_file <- file.path(var_name, paste0("gene_summary_", var_name, ".tsv"))
if (file.exists(gene_summary_file)) {
  gene_summary <- read.table(gene_summary_file, header=TRUE, sep="\t",
                              quote="", comment.char="", stringsAsFactors=FALSE)
  # Every LD-curated candidate gene for this variable, not just the ones
  # flagged In_40000_3_Peak == "Yes" -- that flag comes from a separately
  # precomputed gene_list_<var>_40000_3.txt and doesn't reliably agree
  # with which genes actually have a retained-SNP anchor in THIS script's
  # own live retained-SNP recomputation below (confirmed: several genes
  # the user manually identity-corrected turned out to be flagged "No"
  # here despite genuinely mapping to a retained peak via snp_to_gene).
  # top_snp_per_gene's snp_to_gene-based matching against the live
  # `retained` set is what actually determines inclusion -- a gene with
  # no retained-SNP anchor is silently dropped there regardless.
  genes <- gene_summary
  genes$Gene_Coordinates <- genes$Coordinates
  # A bare Arabidopsis Genome Initiative locus tag (e.g. "At5g09550") isn't
  # a real gene symbol -- it's an anonymous positional ID that happens to
  # satisfy the digit-required acronym-extraction filter upstream, but
  # conveys nothing informative on a plot label, so it's excluded here (and
  # from Orthologs text below) exactly like "uncharacterized protein" is.
  is_agi_locus_id <- function(x) grepl("^At[0-9CM]G[0-9]+$", x, ignore.case = TRUE)

  # Some Acronym values are pure punctuation/whitespace artifacts from
  # upstream parsing (e.g. literal ";" or "; ; ; ; ;") — reject those, but
  # keep legitimate symbols that happen to contain special characters
  # (e.g. "C/VIF1", "EIF(ISO)4E") by requiring at least one alphanumeric.
  has_acronym  <- !is.na(genes$Acronym)   & genes$Acronym   != "" & genes$Acronym   != "NA" &
                   grepl("[A-Za-z0-9]", genes$Acronym) & !is_agi_locus_id(genes$Acronym)

  # Orthologs can be "; "-joined (up to 3 names); some of those are
  # UniProt's generic auto-name for an uncharacterized locus ("AT2G25670
  # protein", "Uncharacterized protein At4g32610") rather than a real gene
  # name — strip those out part-by-part rather than rejecting the whole
  # field, and use it as a label only if something informative remains.
  # Bare classification/process sentences from eggNOG hits that carry no
  # specific gene identity -- e.g. "floral organ abscission" (a GO
  # Biological Process phrase, not a gene name) or "Belongs to the X
  # family" (a UniProt-style classification sentence with nothing named).
  # An exact-match list, not a keyword heuristic, since the full set of
  # eggNOG-sourced Orthologs strings in this dataset is small and has
  # been manually reviewed -- some superficially similar phrasings DO
  # name something real ("Sorting and assembly machinery" = the SAM/TOB
  # complex, "VQ motif" = a real plant protein family) and are
  # deliberately NOT in this list.
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
  # Labels are placed with x fixed at the peak (direction = "y"), so a very
  # wide label can run off the left/right edge of the panel regardless of
  # vertical headroom — keep it to the single best (first-listed) ortholog
  # name rather than joining up to 3, and hard-truncate as a last resort.
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

# ggrepel's own force-directed repulsion only nudges a label off the floor
# when it directly collides with a neighbor — labels with modest x-gaps
# between them get left sitting right at the floor, since that's their
# closest allowed position to the real (much lower) data value. That reads
# as "everything crammed on one line" once there are more than a handful
# of labels. Deterministically assign each label to one of a few height
# tiers (alternating by x-order, so x-adjacent labels land on different
# tiers) via nudge_y, so the vertical room actually gets used; ggrepel
# still does fine local adjustment from each tier's starting point.
assign_tier_nudge <- function(lab_df, y_floor, y_ceiling, y_col, n_tiers = 3) {
  n <- nrow(lab_df)
  if (n == 0) return(numeric(0))
  n_tiers <- min(n_tiers, n)
  tier <- (rank(lab_df$pos, ties.method = "first") - 1) %% n_tiers
  tier_height <- if (n_tiers > 1) (y_ceiling - y_floor) / (n_tiers - 1) else 0
  target_y <- y_floor + tier * tier_height
  target_y - lab_df[[y_col]]
}

# For each gene, find the retained SNP with the highest -log10(p) that
# SnpEff actually assigned it to -- via snp_to_gene_<var_name>.tsv (built
# by build_snp_gene_map.py from outliers_<var>.ann.INFO), NOT a
# coordinate-range overlap against the gene's own Coordinates span. A
# coordinate-only check systematically misses every gene reached only
# through the LD-curated intergenic assignment (curate_intergenic_by_ld.py)
# -- SnpEff's upstream/downstream/intergenic flanking window commonly
# extends several kb outside the gene's own coordinates, which is the
# majority case here (most outlier SNPs are intergenic) and exactly why
# select_peak_labels.R already made this same fix; annotated_manhattan.R
# had inherited the older, narrower coordinate-overlap approach, which is
# why the vast majority of genes flagged In_40000_3_Peak == "Yes" were
# silently ending up unlabeled (41 qualifying genes for BIO15, only 10
# actually plotted, before this fix).
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
  # Match on Gene_ID -- the true key snp_to_gene_<var>.tsv uses -- not
  # Gene_Name, which is only a display label and can differ from Gene_ID
  # (e.g. Gene_ID "IDmodified-gene-78898" displays as Gene_Name
  # "MSTRG.28211"; matching on Gene_Name silently failed for every such
  # gene, since it never appears anywhere in snp_to_gene_<var>.tsv).
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

# Putative inversion region — needed here (not just for the scaffold plot
# below) so it can be excluded from the whole-genome plot's labels.
INV_SCAFFOLD <- "DEINVI_00006"
INV_START <- 43000000
INV_END   <- 50000000
inv_chr <- scaf2chr[[INV_SCAFFOLD]]

# ── Labels: every gene in `candidates` already passed the
# In_40000_3_Peak == "Yes" filter above (applied to `genes` before
# top_snp_per_gene was called), so no further per-peak winner-picking
# happens here — a peak with N qualifying genes gets N labels, not one. ───
label_cols <- c("Chromosome", "Position", "p_final", "Gene_Name", "Label", "has_acronym")
labels <- if (nrow(candidates) > 0) unique(candidates[, label_cols]) else candidates[, label_cols]

# Companion report: which labeled genes are showing an Orthologs-fallback
# label rather than a real acronym — a ready-made worklist for manually
# assigning acronyms to whatever's actually cluttering the plot.
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

# Whole-genome plot excludes labels that fall inside the inversion window —
# those are only shown on the DEINVI_00006 scaffold plot further below.
main_labels <- labels
if (!is.null(inv_chr) && nrow(labels) > 0) {
  in_inversion <- as.character(labels$Chromosome) == inv_chr &
                  labels$Position >= INV_START & labels$Position <= INV_END
  main_labels <- labels[!in_inversion, ]
}

# ── Plot: filtered (40000_3) peaks, gene labels pinned near the top ───────
REPS <- round(length(levels(man$Chromosome)) / 2) + 1
y_top <- max(-log10(man$p_final), na.rm = TRUE)
# geom_text_repel's `ylim` only controls where it PLACES a label — it does
# not feed back into ggplot's automatic axis-range calculation (that's
# still based on each layer's original, un-repelled y aesthetic value), so
# without an explicit ceiling here, labels pushed up near y_top get
# clipped at the top edge of the panel instead of getting more headroom.
# The panel itself gets extra headroom beyond the label ceiling (not just
# the same value) because a label's rendered text glyph has real height
# around its anchor point — clamping the anchor exactly at the panel edge
# still let the tallest labels visually poke past it.
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
  # Labels pinned near the top of the plot (ylim floor, ceiling below the
  # panel's own limit — see y_panel_ceiling note above). Pre-staggered
  # across more height tiers than before (n_tiers=8, not the default 3) —
  # every retained-peak gene is labeled now rather than just an evidence-
  # selected subset, so there are routinely dozens of labels on one plot
  # and 3 tiers left too many crammed on the same line — via nudge_y (see
  # assign_tier_nudge) so the vertical room actually gets used, then
  # ggrepel's own force does fine local adjustment from there (a high
  # box.padding/force pair so it actively separates labels rather than
  # settling for minimal spacing) — direction = "y" keeps them at each
  # peak's own x position, so the connecting segment reads as a straight,
  # thin line up from the peak. All labels are plain black (no GO-term
  # color-coding or legend — dropped per explicit request).
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
# Same per-scaffold Manhattan_40000_3_<scaf>_<var_name> style as
# find_outliers.R, for the scaffold carrying the putative inversion —
# full scaffold, not just the 43-50Mb inversion window — labeled with gene
# names/acronyms (reusing the whole-genome `labels` table computed above —
# the FULL one, not main_labels, so inversion-window genes excluded from
# the whole-genome plot are shown here instead — filtered to this
# scaffold). All labels plain black, same as the whole-genome plot.
# Differences from find_outliers.R's version: y-axis text is doubled
# (theme_bw(base_size=11) default 8.8pt -> 17.6pt here), the Bonferroni
# threshold line is drawn heavier (linewidth 1.5 vs 0.5), and the aspect
# ratio is flatter (8x2.3in vs 7x3.5in). INV_SCAFFOLD/inv_chr were already
# defined above (needed there to exclude inversion labels from the main plot).
if (!is.null(inv_chr)) {
  inv_retained <- retained_df[as.character(retained_df$Chromosome) == inv_chr, ]

  if (nrow(inv_retained) > 0) {
    inv_snps   <- data[as.character(data$Chromosome) == inv_chr, ]
    inv_labels <- labels[as.character(labels$Chromosome) == inv_chr, ]

    FAI <- "/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta.masked.fai"
    fai_df <- read.table(FAI, header = FALSE,
                         col.names = c("scaffold", "length", "offset", "linebases", "linewidth"),
                         stringsAsFactors = FALSE)
    scaf_len <- fai_df$length[fai_df$scaffold == INV_SCAFFOLD]
    x_breaks <- seq(0, scaf_len, by = 10e6)

    # Ceiling matching the label-placement ylim below — see note on
    # y_label_ceiling/y_panel_ceiling in the whole-genome plot above for
    # why these are two different values.
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

print(paste("Done:", var_name))
