#!/usr/bin/env Rscript --vanilla

# Reports EVERY GO category with raw (uncorrected) p < P_NOMINAL_MAX from the
# same goseq gene-length-corrected test GOenrich_goseq.R runs, regardless of
# whether they survive BH/FDR correction. For candidate gene lists this
# small against a GO universe this dense (~9200 testable categories after
# size filtering), FDR correction can leave literally nothing significant
# even when real per-term signal exists (see GOenrich_goseq.R's own log
# output for the min raw p / min FDR for a given run). This is meant for
# reporting nominal hits as hypothesis-generating, explicitly-caveated
# candidates — NOT as a substitute for the FDR-based result.
#
# Because this is now an uncorrected p < 0.05 cut against ~9000 tests (not a
# fixed top-30 list), a lot of what clears that bar is exactly the kind of
# broad, constitutively-used category (translation, proteasome, core
# metabolism, "protein binding") that's virtually guaranteed to contain some
# nominally-enriched genes by chance regardless of the phenotype — see
# build_housekeeping_go_exclusion_list.py for the (deliberately narrow,
# hand-curated) list of such terms filtered out below.
#
# Usage: Rscript top_nominal_GO_hits.R <target_file> [output_dir] [p_max]

if (getRversion() < "4.4.0") `%||%` <- function(x, y) if (is.null(x)) y else x

suppressMessages({
    library(goseq)
    library(dplyr)
    library(GO.db)
    library(rrvgo)
})

MIN_CAT_SIZE <- 3

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript top_nominal_GO_hits.R <target_file> [output_dir] [p_max]")
}
target_file <- args[1]
output_dir  <- if (length(args) >= 2) args[2] else dirname(target_file)
P_NOMINAL_MAX <- if (length(args) >= 3) as.numeric(args[3]) else 0.05

file_label <- sub("^gene_list_|^genes_", "", sub("\\.txt$", "", basename(target_file)))

GFF <- "/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"

# --- Load GO annotation (same combined source as GOenrich_goseq.R) ---
go_annot_raw <- read.table("/private/groups/meyerlab/esteiger/genome/annotated_sequences_combined_geneid_go.tsv",
                           header=TRUE, sep="\t", stringsAsFactors=FALSE)
names(go_annot_raw) <- c("gene", "go")
go_annot_raw$gene <- trimws(go_annot_raw$gene)
go_annot_raw <- unique(go_annot_raw[, c("gene", "go")])

# Remap obsolete GO ids to their current successor where one exists — see
# GSEA_GO.R's comment and build_go_ancestor_map.py for how the remap table
# was built (only unambiguous `replaced_by` chains, not the ambiguous
# `consider:` alternatives).
remap_file <- "/private/groups/meyerlab/esteiger/genome/go_obsolete_remap.tsv"
if (file.exists(remap_file)) {
    remap_tab <- read.table(remap_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    remap_vec <- setNames(remap_tab$new_id, remap_tab$old_id)
    n_remapped <- sum(go_annot_raw$go %in% names(remap_vec))
    go_annot_raw$go <- ifelse(go_annot_raw$go %in% names(remap_vec),
                              remap_vec[go_annot_raw$go], go_annot_raw$go)
    go_annot_raw <- unique(go_annot_raw)
    print(paste0("Remapped ", n_remapped, " obsolete-id gene-GO pairs to their current successor term -> ",
                 nrow(go_annot_raw), " gene-GO pairs after dedup"))
} else {
    print(paste0("WARNING: obsolete-id remap table not found at ", remap_file, " — skipping that step"))
}

# Drop annotations with no plausible plant equivalent (e.g. "embryonic
# skeletal system development") BEFORE propagation and testing, not just
# from the final report — this way excluded terms don't (a) get to seed
# implausible ancestor terms via propagation, or (b) count toward the
# multiple-testing burden. See build_nonplant_go_exclusion_list.py for how
# this list was built and why (BLAST-hit-derived GO transfer isn't
# restricted to plant source organisms). In practice this removes <1% of
# tested categories, so it won't rescue FDR significance on its own — it's
# a data-quality fix, not a power fix.
excl_file <- "/private/groups/meyerlab/esteiger/genome/go_terms_excluded_nonplant.txt"
if (file.exists(excl_file)) {
    excl_ids <- read.table(excl_file, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")[[1]]
    n_before <- nrow(go_annot_raw)
    go_annot_raw <- go_annot_raw[!(go_annot_raw$go %in% excl_ids), ]
    print(paste0("Non-plant-plausible term filter: ", n_before, " -> ", nrow(go_annot_raw), " direct gene-GO pairs"))
} else {
    print(paste0("WARNING: non-plant exclusion list not found at ", excl_file, " — skipping that filter"))
}

# Housekeeping terms are NOT dropped here (unlike the non-plant filter
# above) — they're still real, plausible plant annotations, just not
# informative as a REPORTED enrichment category. Filtering them this early
# would also strip them from propagation/testing entirely, which would
# incorrectly shrink the multiple-testing burden and gene universe for
# every OTHER category too. So this list is loaded now but applied later,
# only at the final per-ontology selection step, after p-values are computed.
housekeeping_file <- "/private/groups/meyerlab/esteiger/genome/go_terms_excluded_housekeeping.txt"
housekeeping_ids <- character(0)
if (file.exists(housekeeping_file)) {
    housekeeping_ids <- read.table(housekeeping_file, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")[[1]]
} else {
    print(paste0("WARNING: housekeeping exclusion list not found at ", housekeeping_file, " — skipping that filter"))
}

# --- Propagate annotations up the GO hierarchy (true path rule) ---
# Ancestors come from a fresh go-basic.obo (build_go_ancestor_map.py), NOT
# GO.db's own GOBPANCESTOR/GOCCANCESTOR/GOMFANCESTOR tables — those were
# found to lag behind the current ontology (missing a recently-added is_a
# edge; see build_nonplant_go_exclusion_list.py for how that was caught).
anc_tab <- read.table("/private/groups/meyerlab/esteiger/genome/go_ancestor_map.tsv",
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
anc_map <- split(anc_tab$ancestor, anc_tab$term)
unique_terms <- unique(go_annot_raw$go)
anc_list <- lapply(unique_terms, function(t) {
    a <- anc_map[[t]]
    if (is.null(a)) return(t)
    a <- a[grepl("^GO:", a)]
    unique(c(t, a))
})
names(anc_list) <- unique_terms

expanded_terms <- anc_list[go_annot_raw$go]
n_anc <- lengths(expanded_terms)
go_annot_prop <- data.frame(
    gene = rep(go_annot_raw$gene, n_anc),
    go   = unlist(expanded_terms, use.names = FALSE),
    stringsAsFactors = FALSE
)
go_annot_prop <- unique(go_annot_prop)

go_annot <- go_annot_prop %>%
  group_by(gene) %>%
  summarise(GO_terms = paste(unique(go), collapse = ",")) %>%
  ungroup()

gene2GO <- strsplit(go_annot$GO_terms, ",")
names(gene2GO) <- go_annot$gene
all_genes <- names(gene2GO)

# --- Read the full GFF once: gene lengths + a transcript_id -> gene_id map ---
gff_all <- read.table(GFF, sep="\t", header=FALSE, quote="", comment.char="#",
                      stringsAsFactors=FALSE,
                      col.names=c("scaffold","source","feature","start","end",
                                  "score","strand","frame","attributes"))

gff_genes <- gff_all[gff_all$feature == "gene", ]
gff_genes$gene_id <- sub(".*gene_id=([^;]+).*", "\\1", gff_genes$attributes)
gene_length <- setNames(gff_genes$end - gff_genes$start + 1, gff_genes$gene_id)
gff_gene_ids <- unique(gff_genes$gene_id)

mrna_rows <- gff_all[gff_all$feature == "mRNA", ]
tx_ids <- sub(".*ID=([^;]+).*", "\\1", mrna_rows$attributes)
gene_ids_col <- sub(".*Parent=([^;]+).*", "\\1", mrna_rows$attributes)
has_both <- grepl("ID=", mrna_rows$attributes) & grepl("Parent=", mrna_rows$attributes)
tx2gene_gff <- setNames(gene_ids_col[has_both], tx_ids[has_both])
tx2gene_gff <- tx2gene_gff[!duplicated(names(tx2gene_gff))]

# --- Load target gene list ---
target_genes <- unique(trimws(scan(target_file, what="")))
target_genes <- sub("\\.$", "", target_genes)

resolve_one <- function(g) {
    cand <- c(g, paste0(g, ".t1"), paste0(g, ".1"))
    hit <- cand[cand %in% names(tx2gene_gff)]
    if (length(hit) > 0) return(unname(tx2gene_gff[hit[1]]))
    pref <- grep(paste0("^", gsub("\\.", "\\\\.", g), "\\."), names(tx2gene_gff), value = TRUE)
    if (length(pref) > 0) return(unname(tx2gene_gff[pref[1]]))
    NA_character_
}
not_a_gene_id <- !(target_genes %in% gff_gene_ids)
resolvable <- vapply(target_genes[not_a_gene_id], resolve_one, character(1))
n_resolved <- sum(!is.na(resolvable))
if (n_resolved > 0) {
    target_genes[not_a_gene_id][!is.na(resolvable)] <- resolvable[!is.na(resolvable)]
    target_genes <- unique(target_genes)
}

targets <- length(target_genes)
conserved <- length(target_genes[target_genes %in% all_genes])
print(paste0("Number of targets was ", targets, ", but number found in gene universe was ", conserved))
target_genes <- target_genes[target_genes %in% all_genes]

# --- Align universe to genes with known length ---
common_genes <- intersect(all_genes, names(gene_length))
gene.vector <- as.integer(common_genes %in% target_genes)
names(gene.vector) <- common_genes
bias.data <- gene_length[common_genes]

gene2cat <- gene2GO[common_genes]
cat_counts <- table(unlist(gene2cat))
keep_cats <- names(cat_counts)[cat_counts >= MIN_CAT_SIZE]
gene2cat <- lapply(gene2cat, function(x) x[x %in% keep_cats])
print(paste0("Category size filter (>=", MIN_CAT_SIZE, " genes): ",
             length(cat_counts), " categories -> ", length(keep_cats), " kept"))

pwf <- nullp(gene.vector, bias.data = bias.data, plot.fit = FALSE)

GO.wall <- goseq(pwf, gene2cat = gene2cat, test.cats = c("GO:BP", "GO:MF", "GO:CC"))
GO.wall$FDR <- p.adjust(GO.wall$over_represented_pvalue, method = "BH")
GO.wall <- GO.wall[order(GO.wall$over_represented_pvalue), ]

# --- Select nominal hits (raw p < P_NOMINAL_MAX), drop housekeeping terms,
# then reduce redundancy via rrvgo, per ontology ---
# Pooling and the final report use the SAME p < P_NOMINAL_MAX cut (unlike
# the old top-30-by-rank version, which had to pool well beyond n_top so
# reduction wouldn't crowd out terms that just missed the rank cutoff — that
# headroom isn't needed once selection is threshold-based instead of
# rank-based, since the pool already IS exactly the set eligible for the
# final report). Housekeeping terms are dropped from the pool before rrvgo
# runs, so they can't crowd out real signal there either, not just from the
# final table. Redundancy reduction itself has the same rationale as
# GOenrich_goseq.R's rrvgo step for FDR-significant terms: without it,
# chains of nested parent/child terms (all driven by the same underlying
# gene(s) via propagation) would look like several independent hits when
# they're substantially one signal. Scored by raw p here (not FDR, which is
# 1 almost everywhere for nominal-only hits — see module docstring).
# method="Wang" is graph-topology-only (no annotation-frequency/Information
# Content data needed), appropriate since org.At.tair.db is only a stand-in
# orgdb for this species.
reduced_list <- list()
n_pooled <- 0
n_housekeeping_dropped <- 0
for (GROUP in c("BP", "MF", "CC")) {
    sub_res <- GO.wall[!is.na(GO.wall$ontology) & GO.wall$ontology == GROUP, ]
    sub_res <- sub_res[sub_res$over_represented_pvalue < P_NOMINAL_MAX, ]
    n_housekeeping_dropped <- n_housekeeping_dropped + sum(sub_res$category %in% housekeeping_ids)
    sub_res <- sub_res[!(sub_res$category %in% housekeeping_ids), ]
    n_pooled <- n_pooled + nrow(sub_res)
    if (nrow(sub_res) >= 3) {
        simMatrix <- calculateSimMatrix(sub_res$category, orgdb = "org.At.tair.db",
                                        ont = GROUP, method = "Wang")
        reducedTerms <- reduceSimMatrix(simMatrix,
                                        scores = setNames(-log10(sub_res$over_represented_pvalue), sub_res$category),
                                        threshold = 0.7, orgdb = "org.At.tair.db")
        reduced_GO <- reducedTerms[reducedTerms$parent == reducedTerms$go, ]
        sub_res <- sub_res[sub_res$category %in% reduced_GO$go, ]
    }
    reduced_list[[GROUP]] <- sub_res
}
GO.wall.reduced <- do.call(rbind, reduced_list)
GO.wall.reduced <- GO.wall.reduced[order(GO.wall.reduced$over_represented_pvalue), ]
print(paste0("Nominal hits (raw p < ", P_NOMINAL_MAX, "), housekeeping terms dropped: ",
             n_housekeeping_dropped, " | after rrvgo redundancy reduction: ",
             n_pooled, " -> ", nrow(GO.wall.reduced), " categories"))

gene_by_cat <- split(rep(names(gene2cat), lengths(gene2cat)), unlist(gene2cat))
top <- GO.wall.reduced
top$Genes <- sapply(top$category, function(go_id) {
    genes <- gene_by_cat[[go_id]]
    cand_genes <- genes[genes %in% target_genes]
    paste(cand_genes, collapse = ", ")
})

outfile <- file.path(output_dir, paste0("GO_top_nominal_", file_label, ".tsv"))
write.table(top, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
print(paste0("Wrote ", nrow(top), " nominal hits (raw p < ", P_NOMINAL_MAX,
             ", NOT FDR-significant, housekeeping terms excluded) to ", outfile))
print(paste0("Min raw p-value: ", format(min(GO.wall$over_represented_pvalue, na.rm=TRUE), digits=3),
             "  |  Min FDR: ", format(min(GO.wall$FDR, na.rm=TRUE), digits=3),
             "  |  Categories tested: ", nrow(GO.wall)))
