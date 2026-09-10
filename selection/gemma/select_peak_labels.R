#!/usr/bin/env Rscript
#
# Reads the kinship-corrected mixed-model GEMMA output (output/pheno<N>_lmm_
# results.assoc.txt), NOT the older un-corrected output/pheno<N>_results.assoc.txt
# (deleted — it let a large chunk of the DEINVI_00006 44-49Mb inversion region's
# structure-confounded signal through: 373 "Bonferroni-significant" SNPs there
# under the old model collapsed to 3 once kinship was modeled properly).
#
# For each independent peak in the 40000_3 retained-SNP set (same Bonferroni
# + 40kb-window retention + 40kb peak-merge logic as annotated_manhattan.R),
# select a single gene to label it with, combining two factors via a
# rank-sum (lower total rank wins; GEMMA p-value and mutation impact are
# NOT used to choose which gene is selected, only carried through as
# informational columns):
#   (A) GSEA GO-term enrichment — is the gene in the leading edge of any
#       significantly enriched Biological Process (BP) term (pval < 0.005)
#       for this variable?
#   (C) Functional tie to adaptation — Putative_Adaptive_Role/Citation from
#       gene_summary_all.tsv (real citation > role but no citation > none).
#
# A peak's candidate genes are found via <var_name>/snp_to_gene_<var_name>.tsv
# (built by build_snp_gene_map.py from outliers_<var_name>.ann.INFO using the
# same 3-tier ID resolver validated for the Mutation_Types/Gene_Coordinates
# fix), i.e. each retained SNP's ACTUAL SnpEff-assigned gene(s) — not
# coordinate-range overlap, which misses genes SnpEff correctly assigned via
# its upstream/downstream/intergenic flanking window (commonly several kb
# outside the gene's own coordinate span, and the majority case here since
# ~78% of outlier SNPs are intergenic).
#
# Writes <var_name>/for_manhattan_label_<var_name>.tsv with only the winning
# gene per peak (one row per peak that has >=1 labelable candidate; peaks
# with none are simply omitted).
#
# Only genes with a usable label are considered at all (same Acronym ->
# cleaned-Ortholog -> exclude logic as annotated_manhattan.R) — there'd be
# nothing to write on the plot otherwise.

args <- commandArgs(trailingOnly=TRUE)
pheno_num <- args[1]
var_name  <- args[2]

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs")
library(data.table)
library(dplyr)

setwd("/private/groups/meyerlab/esteiger/selection/gemma/")
dir.create(var_name, showWarnings=FALSE, recursive=TRUE)

dic <- read.table("scaffold_snp_map.txt")
dic <- unique(dic[,c("V1","V3")])
scaf2chr <- setNames(as.character(dic$V1), dic$V3)
chr2scaf <- setNames(dic$V3, as.character(dic$V1))

# ── Recompute retained (40000_3) SNPs — identical to annotated_manhattan.R ──
gemma <- read.table(paste0("output/pheno", pheno_num, "_lmm_results.assoc.txt"), header=T)
names(gemma)[1:3] <- c("Chromosome", "SNP", "Position")
print(paste("No rows in input:", nrow(gemma)))

chisq_obs <- qchisq(gemma$p_lrt, df=1, lower.tail=FALSE)
lambda <- median(chisq_obs, na.rm=TRUE) / qchisq(0.5, df=1)
if (lambda > 1.05) {
  gemma$p_final <- pchisq(chisq_obs / lambda, df=1, lower.tail=FALSE)
} else {
  gemma$p_final <- gemma$p_lrt
}
thresh_bonf <- 0.05 / nrow(gemma)
pks <- as.data.table(gemma[gemma$p_final < thresh_bonf, ])
print(paste("Bonferroni hits:", nrow(pks)))

if (nrow(pks) == 0) {
  print(paste("No Bonferroni outliers for", var_name, "— nothing to select"))
  quit(save="no", status=0)
}

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
setnames(win_kept, c("start", "end"), c("wstart", "wend"), skip_absent = TRUE)

pks[, pos2 := Position]
setkey(pks, Chromosome, Position, pos2)
setkey(win_kept, Chromosome, wstart, wend)
retained <- unique(foverlaps(pks, win_kept, by.x = c("Chromosome", "Position", "pos2"),
                      by.y = c("Chromosome", "wstart", "wend"), nomatch = 0L))
retained[, c("wstart", "wend", "N", "pos2") := NULL]
retained <- unique(retained)
print(paste("Retained SNPs (multi-SNP peaks):", nrow(retained)))

if (nrow(retained) == 0) {
  print(paste("No 40000_3 peaks for", var_name, "— nothing to select"))
  quit(save="no", status=0)
}
retained_df <- as.data.frame(retained)
retained_df$Chromosome <- as.character(retained_df$Chromosome)

# ── Cluster retained SNPs into independent peaks (40kb merge gap) ──────────
PEAK_MERGE_GAP <- 40000
retained_df <- retained_df[order(retained_df$Chromosome, retained_df$Position), ]
peak_id <- integer(nrow(retained_df)); cur <- 0; prev_chr <- NA; prev_pos <- -Inf
for (i in seq_len(nrow(retained_df))) {
  if (!identical(retained_df$Chromosome[i], prev_chr) || (retained_df$Position[i] - prev_pos) > PEAK_MERGE_GAP) cur <- cur + 1
  peak_id[i] <- cur
  prev_chr <- retained_df$Chromosome[i]; prev_pos <- retained_df$Position[i]
}
retained_df$peak_id <- peak_id

peaks <- retained_df %>%
  group_by(peak_id, Chromosome) %>%
  summarise(Peak_Start = min(Position), Peak_End = max(Position), Peak_N_SNPs = n(), .groups = "drop")
peaks$Scaffold <- chr2scaf[peaks$Chromosome]
print(paste("Independent peaks:", nrow(peaks)))

# ── Gene name/label lookup — same Acronym -> cleaned-Ortholog -> exclude
# logic as annotated_manhattan.R ────────────────────────────────────────────
gene_summary <- read.table("gene_summary_all.tsv", header=TRUE, sep="\t",
                            quote="", comment.char="", stringsAsFactors=FALSE)
genes <- gene_summary[gene_summary$Variable == var_name, ]
# Some Acronym values are pure punctuation/whitespace artifacts from
# upstream parsing (e.g. literal ";" or "; ; ; ; ;") — reject those, but
# keep legitimate symbols that happen to contain special characters
# (e.g. "C/VIF1", "EIF(ISO)4E") by requiring at least one alphanumeric.
has_acronym <- !is.na(genes$Acronym) & genes$Acronym != "" & genes$Acronym != "NA" &
                grepl("[A-Za-z0-9]", genes$Acronym)
is_generic_ortholog <- function(x) {
  grepl("^AT[0-9CM]G[0-9]+(\\s+protein)?$", x, ignore.case = TRUE) |
  grepl("uncharacterized protein", x, ignore.case = TRUE)
}
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
genes$Label <- ifelse(has_acronym, genes$Acronym, genes$Orthologs_clean)
genes <- genes[!is.na(genes$Label), ]
print(paste("Genes for", var_name, "with a usable label:", nrow(genes)))

# ── SNP -> gene mapping (built by build_snp_gene_map.py from
# outliers_<var>.ann.INFO, using the same 3-tier ID resolver validated for
# the Mutation_Types/Gene_Coordinates fix) — see header note for why this
# replaces coordinate-range overlap ────────────────────────────────────────
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

retained_df$Scaffold <- chr2scaf[retained_df$Chromosome]

# gene -> data.frame(peak_id, p_final) for every retained SNP resolved to it
gene_hits <- new.env()
for (i in seq_len(nrow(retained_df))) {
  gids <- genes_for_snp(retained_df$Scaffold[i], retained_df$Position[i])
  if (length(gids) == 0) next
  row <- data.frame(peak_id = retained_df$peak_id[i], p_final = retained_df$p_final[i])
  for (gid in gids) {
    cur <- get0(gid, envir = gene_hits, ifnotfound = NULL)
    assign(gid, rbind(cur, row), envir = gene_hits)
  }
}

# ── GSEA leading-edge lookup: gene -> significant terms it drives ─────────
# GSEA_GO_results_<var>.tsv is already filtered to nominal p<0.05, but genes
# are only credited with a GO-term hit (factor A) for Biological Process
# (BP) terms stricter than that, at pval < 0.005, per user request.
GSEA_PVAL_CUTOFF <- 0.005
gsea_file <- file.path(var_name, paste0("GSEA_GO_results_", var_name, ".tsv"))
gene_go_terms <- new.env()
if (file.exists(gsea_file)) {
  gsea <- read.table(gsea_file, header=TRUE, sep="\t", quote="", comment.char="",
                      stringsAsFactors=FALSE)
  gsea <- gsea[!is.na(gsea$pval) & gsea$pval < GSEA_PVAL_CUTOFF &
               !is.na(gsea$ontology) & gsea$ontology == "BP", ]
  for (i in seq_len(nrow(gsea))) {
    lead <- trimws(strsplit(gsea$leadingEdge[i], ",")[[1]])
    for (g in lead) {
      cur <- get0(g, envir = gene_go_terms, ifnotfound = character(0))
      assign(g, c(cur, gsea$term[i]), envir = gene_go_terms)
    }
  }
  print(paste("Significant BP GSEA terms loaded (p <", GSEA_PVAL_CUTOFF, "):", nrow(gsea)))
} else {
  print(paste("No GSEA results file found for", var_name, "— factor A will be empty"))
}
gsea_terms_for <- function(gid) {
  terms <- get0(gid, envir = gene_go_terms, ifnotfound = character(0))
  if (length(terms) == 0) return(NA_character_)
  paste(unique(terms), collapse = "; ")
}

# ── Impact severity ranking ─────────────────────────────────────────────
impact_order <- c(HIGH = 0, MODERATE = 1, LOW = 2, MODIFIER = 3)
best_impact <- function(impact_str) {
  if (is.na(impact_str) || impact_str == "NA" || impact_str == "") return(NA_character_)
  levels <- trimws(strsplit(impact_str, ";")[[1]])
  levels <- levels[levels %in% names(impact_order)]
  if (length(levels) == 0) return(NA_character_)
  levels[which.min(impact_order[levels])]
}

# ── Build candidate table: for each labelable gene, use gene_hits (built
# above from the SNP -> gene map) to find which peak(s) it was actually
# SnpEff-assigned to, and its best p-value within each ────────────────────
candidates <- list()
for (i in seq_len(nrow(genes))) {
  gid <- genes$Gene_Name[i]
  hits <- get0(gid, envir = gene_hits, ifnotfound = NULL)
  if (is.null(hits)) next
  has_role <- !is.na(genes$Putative_Adaptive_Role[i]) & genes$Putative_Adaptive_Role[i] != "NA"
  has_citation <- !is.na(genes$Citation[i]) & genes$Citation[i] != "NA" &
                   !startsWith(genes$Citation[i], "NOT FOUND")
  adaptive_status <- if (has_role) {
    if (has_citation) "Role+Citation" else "Role only"
  } else "None"
  imp <- best_impact(genes$Impact_Levels[i])
  for (pid in unique(hits$peak_id)) {
    p_in_peak <- min(hits$p_final[hits$peak_id == pid])
    candidates[[length(candidates) + 1]] <- data.frame(
      peak_id = pid,
      Gene_ID = gid,
      Label = genes$Label[i],
      Best_P_Value_in_Peak = p_in_peak,
      GSEA_GO_Terms = gsea_terms_for(gid),
      Adaptive_Role_Status = adaptive_status,
      Putative_Adaptive_Role = genes$Putative_Adaptive_Role[i],
      Best_Impact_Level = ifelse(is.na(imp), "NA", imp),
      Mutation_Types = genes$Mutation_Types[i],
      stringsAsFactors = FALSE
    )
  }
}

out_cols <- c("Scaffold","Position_Range","Peak_N_SNPs",
              "Label","Best_P_Value_in_Peak",
              "GSEA_Hit","GSEA_GO_Terms","Adaptive_Role_Status","Putative_Adaptive_Role",
              "Best_Impact_Level","Mutation_Types")

if (length(candidates) == 0) {
  print(paste("No labelable candidate genes found for any peak in", var_name))
  cand_df <- data.frame()
} else {
  cand_df <- do.call(rbind, candidates)

  # ── Rank-sum selection within each peak — GO-term association (BP,
  # p<0.005) and adaptive-role support only; GEMMA p-value and mutation
  # impact are not used to pick the winner (informational columns only).
  # Ties are broken alphabetically by Label for a deterministic pick.
  adaptive_badness <- c("Role+Citation" = 0, "Role only" = 1, "None" = 2)
  cand_df$GSEA_Hit   <- !is.na(cand_df$GSEA_GO_Terms)
  cand_df$A_badness  <- ifelse(cand_df$GSEA_Hit, 0, 1)
  cand_df$C_badness  <- adaptive_badness[cand_df$Adaptive_Role_Status]

  cand_df <- cand_df %>%
    group_by(peak_id) %>%
    mutate(
      Rank_A   = rank(A_badness, ties.method = "min"),
      Rank_C   = rank(C_badness, ties.method = "min"),
      Rank_Sum = Rank_A + Rank_C
    ) %>%
    arrange(peak_id, Rank_Sum, Label) %>%
    mutate(Selected = row_number() == 1) %>%
    ungroup() %>%
    as.data.frame()

  peaks_with_cand <- unique(cand_df$peak_id)

  # Keep only the winning gene per peak (Selected == TRUE); the full
  # per-peak candidate breakdown (all genes considered, individual factor
  # ranks) is available by removing this filter if the audit trail is
  # needed again.
  cand_df <- cand_df[cand_df$Selected == TRUE, ]
  cand_df <- merge(cand_df, peaks[, c("peak_id","Scaffold","Peak_Start","Peak_End","Peak_N_SNPs")], by = "peak_id")

  # If the same gene wins more than one peak (common for adjacent 40kb
  # windows both falling in that gene's flanking SnpEff-assignment range),
  # collapse to a single row spanning all of those peaks' coordinates
  # rather than repeating the gene label on the plot.
  cand_df <- cand_df %>%
    group_by(Gene_ID) %>%
    summarise(
      Scaffold = dplyr::first(Scaffold),
      Peak_Start = min(Peak_Start),
      Peak_End = max(Peak_End),
      Peak_N_SNPs = sum(Peak_N_SNPs),
      Label = dplyr::first(Label),
      Best_P_Value_in_Peak = min(Best_P_Value_in_Peak),
      GSEA_Hit = any(GSEA_Hit),
      GSEA_GO_Terms = dplyr::first(GSEA_GO_Terms),
      Adaptive_Role_Status = dplyr::first(Adaptive_Role_Status),
      Putative_Adaptive_Role = dplyr::first(Putative_Adaptive_Role),
      Best_Impact_Level = dplyr::first(Best_Impact_Level),
      Mutation_Types = dplyr::first(Mutation_Types),
      .groups = "drop"
    ) %>%
    as.data.frame()

  cand_df$Position_Range <- paste0(cand_df$Peak_Start, "-", cand_df$Peak_End)
  cand_df <- cand_df[order(cand_df$Scaffold, cand_df$Peak_Start), ]
  cand_df <- cand_df[, out_cols]
}

if (!exists("peaks_with_cand")) peaks_with_cand <- integer(0)

out_file <- file.path(var_name, paste0("for_manhattan_label_", var_name, ".tsv"))
write.table(cand_df, out_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

print(paste("Peaks with >=1 labelable candidate:", length(peaks_with_cand), "/", nrow(peaks)))
print(paste("Wrote", nrow(cand_df), "selected rows to", out_file))
print(paste("Done:", var_name))
