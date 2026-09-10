#!/usr/bin/env Rscript --vanilla

# GSEA-style GO enrichment: ranks EVERY gene by a continuous statistic
# derived from its SNPs' GEMMA p-values (build_gene_gwas_scores.py's
# ACAT-combined per-gene p-value), and tests for enrichment along that
# ranking with fgsea — as opposed to GOenrich_goseq.R's binary candidate/
# non-candidate split built by thresholding SNPs at Bonferroni significance.
# That binary approach caps candidate lists at a few hundred genes against
# ~9000 tested categories, which is why FDR correction there finds nothing
# significant for any variable even when real per-term signal may exist.
# Using the continuous p-values for every gene (not just the ones that
# cleared a hard cutoff) gives fgsea much more to work with.
#
# Usage: Rscript GSEA_GO.R <var_name> [output_dir]
#   e.g. Rscript GSEA_GO.R fog fog
# Expects <var_name>_gene_gwas_scores.tsv (build_gene_gwas_scores.py output)
# in output_dir (default: same directory this is run from, matching
# build_gene_gwas_scores.py's own default output location).

suppressMessages({
    library(dplyr)
    library(fgsea)
})

P_MAX <- 0.05

MIN_CAT_SIZE <- 3

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript GSEA_GO.R <var_name> [output_dir] [go_annot_file] [output_tag]")
}
var_name    <- args[1]
output_dir  <- if (length(args) >= 2) args[2] else var_name
# go_annot_file lets an alternate gene-GO annotation source (e.g. an
# orthology-only one) be swapped in without duplicating this whole script.
# output_tag, if given, is appended to every output filename (e.g.
# "orthofinder_hannuus" -> GSEA_GO_results_<var>_orthofinder_hannuus.tsv) so
# an alternate-source run never collides with the default combined-source
# results for the same variable.
go_annot_file <- if (length(args) >= 3) args[3] else
    "/private/groups/meyerlab/esteiger/genome/annotated_sequences_combined_geneid_go.tsv"
output_tag  <- if (length(args) >= 4) args[4] else ""
tag_suffix  <- if (output_tag != "") paste0("_", output_tag) else ""

scores_file <- file.path(output_dir, paste0(var_name, "_gene_gwas_scores.tsv"))

# --- Load GO annotation (same propagation as GOenrich_goseq.R regardless of source) ---
go_annot_raw <- read.table(go_annot_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
names(go_annot_raw) <- c("gene", "go")
go_annot_raw$gene <- trimws(go_annot_raw$gene)
go_annot_raw <- unique(go_annot_raw[, c("gene", "go")])
print(paste0("Loaded ", nrow(go_annot_raw), " gene-GO pairs across ",
             length(unique(go_annot_raw$gene)), " genes"))

# Remap obsolete GO ids to their current successor where one exists — the
# combined annotation source was built against whatever GO release eggNOG/
# NCBI/TAIR used at the time, and ~12% of its distinct GO ids have since
# been made obsolete (checked directly: 1,488 of 12,500). An obsolete id
# can't be propagated to ancestors (no entry in go_ancestor_map.tsv) and
# has no display name, so without this it's tested as an orphaned,
# under-connected category. See build_go_ancestor_map.py for how the remap
# table was built (only unambiguous `replaced_by` chains, resolved
# transitively to a live term — NOT the ambiguous `consider:` alternatives).
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

# Drop non-plant-plausible terms BEFORE propagation/testing — union of two
# independently-built exclusion lists: build_nonplant_go_exclusion_list.py's
# hand-picked anchor-term heuristic (broad coverage, but only as good as the
# anchors chosen), and build_taxon_constraint_go_exclusion_list.py's set,
# anchored on the GO Consortium's own curated never_in_taxon-Viridiplantae
# assertions (narrower, but real curation rather than a guess — and catches
# real gaps the heuristic missed entirely, e.g. the whole "centrosome"
# family, which isn't reachable from any of the heuristic's organ/tissue-
# level anchors since centrosomes aren't an organ-system concept at all).
excl_file       <- "/private/groups/meyerlab/esteiger/genome/go_terms_excluded_nonplant.txt"
excl_file_taxon <- "/private/groups/meyerlab/esteiger/genome/go_terms_excluded_taxon_constraint.txt"
excl_ids <- character(0)
if (file.exists(excl_file)) {
    excl_ids <- c(excl_ids, read.table(excl_file, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")[[1]])
} else {
    print(paste0("WARNING: non-plant exclusion list not found at ", excl_file, " — skipping that filter"))
}
if (file.exists(excl_file_taxon)) {
    excl_ids <- c(excl_ids, read.table(excl_file_taxon, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")[[1]])
} else {
    print(paste0("WARNING: taxon-constraint exclusion list not found at ", excl_file_taxon, " — skipping that filter"))
}
excl_ids <- unique(excl_ids)
if (length(excl_ids) > 0) {
    n_before <- nrow(go_annot_raw)
    go_annot_raw <- go_annot_raw[!(go_annot_raw$go %in% excl_ids), ]
    print(paste0("Non-plant-plausible term filter (heuristic + taxon-constraint union, ",
                 length(excl_ids), " ids): ", n_before, " -> ", nrow(go_annot_raw), " direct gene-GO pairs"))
}

# Also drop generic "housekeeping" terms (translation, proteasome, core
# metabolism, "protein binding", ontology roots, etc.) BEFORE testing —
# same list used to clean up top_nominal_GO_hits.R's report
# (build_housekeeping_go_exclusion_list.py), but applied at a different
# stage here deliberately. top_nominal_GO_hits.R has no real multiple-
# testing budget to protect (it just reports raw p < 0.05, uncorrected), so
# it filters housekeeping terms from the report only, after testing, to
# keep an honest "categories tested" count. GSEA_GO.R's actual significance
# criterion IS a multiple-testing correction (fgsea's BH-adjusted padj
# across every tested pathway) — leaving housekeeping/near-universal terms
# (e.g. GO:0005575 "cellular_component", an ontology root) in the tested
# pool spends part of that correction's budget on categories we already
# know a priori aren't informative, diluting power for the ones that are.
# This is legitimate pre-filtering, not data-dependent cherry-picking: the
# exclusion list is fixed and built independently of any particular GEMMA
# run's results.
housekeeping_file <- "/private/groups/meyerlab/esteiger/genome/go_terms_excluded_housekeeping.txt"
if (file.exists(housekeeping_file)) {
    housekeeping_ids <- read.table(housekeeping_file, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")[[1]]
    n_before <- nrow(go_annot_raw)
    go_annot_raw <- go_annot_raw[!(go_annot_raw$go %in% housekeeping_ids), ]
    print(paste0("Housekeeping term filter: ", n_before, " -> ", nrow(go_annot_raw), " direct gene-GO pairs"))
} else {
    print(paste0("WARNING: housekeeping exclusion list not found at ", housekeeping_file, " — skipping that filter"))
}

# --- Propagate annotations up the GO hierarchy (true path rule) ---
# Same fresh-go-basic.obo ancestor map as GOenrich_goseq.R/top_nominal_GO_hits.R
# (NOT GO.db's own ancestor tables — see build_go_ancestor_map.py for why).
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
print(paste0("GO propagation: ", nrow(go_annot_raw), " direct gene-term pairs -> ",
             nrow(go_annot_prop), " after adding ancestors"))

# Housekeeping terms can come BACK after propagation even though they were
# filtered from direct annotations above: propagation walks every surviving
# direct term up to its ontology root, and ontology-root/near-root terms
# (e.g. GO:0005575 "cellular_component") are ancestors of huge swaths of
# otherwise-legitimate, specific terms. Filtering direct annotations alone
# doesn't stop them from being reconstructed as someone else's ancestor —
# confirmed empirically: GO:0005575, GO:0110165, and 6 other housekeeping
# IDs still showed up in the tested pathway set before this second pass was
# added. Re-applying the same filter after propagation catches that.
if (file.exists(housekeeping_file)) {
    n_before <- nrow(go_annot_prop)
    go_annot_prop <- go_annot_prop[!(go_annot_prop$go %in% housekeeping_ids), ]
    print(paste0("Housekeeping term filter (post-propagation): ", n_before, " -> ", nrow(go_annot_prop), " gene-term pairs"))
}
# Same defensive re-pass for non-plant-plausible terms — belt-and-suspenders
# alongside the housekeeping one above. Less likely to matter in practice
# (the non-plant anchors are narrow, specific branches rather than
# near-universal ontology roots, so they're less prone to being someone
# else's ancestor), but cheap to guard against regardless.
if (length(excl_ids) > 0) {
    n_before <- nrow(go_annot_prop)
    go_annot_prop <- go_annot_prop[!(go_annot_prop$go %in% excl_ids), ]
    print(paste0("Non-plant-plausible term filter (post-propagation): ", n_before, " -> ", nrow(go_annot_prop), " gene-term pairs"))
}

go_annot <- go_annot_prop %>%
  group_by(gene) %>%
  summarise(GO_terms = paste(unique(go), collapse = ",")) %>%
  ungroup()

gene2GO <- strsplit(go_annot$GO_terms, ",")
names(gene2GO) <- go_annot$gene

# --- Category size filter (>= MIN_CAT_SIZE genes with any GO annotation) ---
# Same nodeSize-style pruning as GOenrich_goseq.R, but sized against the full
# annotated gene set here (not a length-matched universe — fgsea doesn't need
# gene-length bias correction the way goseq's Wallenius test does, since
# every gene is ranked, not split into a binary in/out test against a
# length-matched background).
cat_counts <- table(unlist(gene2GO))
keep_cats <- names(cat_counts)[cat_counts >= MIN_CAT_SIZE]
gene2cat <- lapply(gene2GO, function(x) x[x %in% keep_cats])
pathways <- split(rep(names(gene2cat), lengths(gene2cat)), unlist(gene2cat))
print(paste0("Category size filter (>=", MIN_CAT_SIZE, " genes): ",
             length(cat_counts), " categories -> ", length(pathways), " kept"))

# --- Load gene-level ACAT p-values and build the ranked stat vector ---
scores <- read.table(scores_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(paste0("Loaded ACAT gene scores for ", nrow(scores), " genes from ", scores_file))

scores <- scores[scores$gene_id %in% names(gene2cat), ]
print(paste0(nrow(scores), " of those genes have GO annotation — the rest can't contribute to any category"))

# Rank by snp_count_adjusted_score (build_gene_gwas_scores.py): -log10(ACAT p)
# with a per-gene LOWESS correction for its (log) SNP count removed. ACAT's
# null distribution is theoretically invariant to SNP count, but checked
# empirically across all 3 variables this pipeline currently runs, a real
# Spearman correlation between n_snps and -log10(ACAT p) showed up for 2 of
# 3 (fog/BIO15, both resolved to ~0 after this correction; BIO18 had none to
# begin with, unaffected) — likely phenotype-specific residual structure
# rather than a flaw in ACAT itself, but worth removing regardless. The
# adjusted score is a signed residual (can be negative), so it's rank-
# transformed to a plain 1..N score before use — fgsea's scoreType="pos"
# below needs non-negative values, and rank() preserves the exact ordering
# that's all a preranked test actually uses.
stats <- rank(scores$snp_count_adjusted_score)
names(stats) <- scores$gene_id
stats <- sort(stats, decreasing = TRUE)

# --- Run fgsea ---
# scoreType = "pos": the ranking stat (-log10 p) is unsigned/single-tailed —
# every gene's value is >= 0, with higher meaning more associated. There's no
# "down" direction the way there is for a signed statistic like differential-
# expression log fold-change, which is what fgsea's default scoreType="std"
# (two-tailed) assumes. Using the default here produced miscalibrated
# results (a "values are never negative" warning, plus 3 pathways where
# p-values couldn't be computed at all due to "unbalanced" pos/neg stats).
set.seed(1)
fgsea_res <- fgsea(pathways = pathways, stats = stats, minSize = MIN_CAT_SIZE,
                    scoreType = "pos", eps = 0)
fgsea_res <- as.data.frame(fgsea_res)  # fgsea returns a data.table; force plain
                                        # data.frame semantics for the column
                                        # selection/reordering below (data.table's
                                        # `dt[, char_vector]` is NOT a column
                                        # subset the way `df[, char_vector]` is)
fgsea_res <- fgsea_res[order(fgsea_res$pval), ]
fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, paste, collapse = ", ")
n_tested <- nrow(fgsea_res)

# --- Human-readable term names + ontology (BP/MF/CC) ---
# From go_term_names.tsv (build_go_ancestor_map.py's fresh-obo parse), NOT
# GO.db's own TERM table — GO.db was found to be missing names for ~40-56
# terms per run here (e.g. GO:0170060 "microtubule destabilizing activity"),
# the same underlying staleness gap already documented for GO.db's ancestor
# tables elsewhere in this pipeline (build_go_ancestor_map.py,
# build_nonplant_go_exclusion_list.py). A pathway with no match here is
# almost always a GO ID that's since been made obsolete in the current
# ontology (still present in the combined annotation source, which was
# built against whatever GO release eggNOG/NCBI/TAIR used at the time) —
# term/ontology both come back NA for those, same root cause either way.
term_name_tab <- read.table("/private/groups/meyerlab/esteiger/genome/go_term_names.tsv",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "",
                            comment.char = "")  # some term names literally contain "#"
                                                # (e.g. GO:1904067 "ascr#2 binding"), which
                                                # read.table's default comment.char="#"
                                                # would otherwise truncate mid-line
match_idx <- match(fgsea_res$pathway, term_name_tab$id)
fgsea_res$term <- term_name_tab$name[match_idx]
fgsea_res$ontology <- term_name_tab$ontology[match_idx]
fgsea_res <- fgsea_res[, c("pathway", "term", "ontology", setdiff(names(fgsea_res), c("pathway", "term", "ontology")))]

# --- Report only nominal hits (raw p < P_MAX), same threshold convention as
# top_nominal_GO_hits.R --- (housekeeping terms are already excluded from
# `pathways` itself, before fgsea ran — see the filter near the top of this
# script — so there's nothing further to drop here.)
fgsea_res <- fgsea_res[!is.na(fgsea_res$pval) & fgsea_res$pval < P_MAX, ]

# --- Annotate leadingEdge genes with a gene acronym (preferred) or,
# failing that, ortholog name, per gene_summary_all.tsv (a combined summary
# across variables' candidate lists, compile_gene_summary.py-style) ---
# NOTE: this file covers each variable's own (threshold-based) candidate
# gene set, not the full ~24,000-gene GSEA universe — so most leadingEdge
# genes here will still have no match and stay as a bare gene id.
summary_file <- "/private/groups/meyerlab/esteiger/selection/gemma/gene_summary_all.tsv"
gene_ortho_label <- character(0)
if (file.exists(summary_file)) {
    gs <- read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                     quote = "", comment.char = "")
    # A gene can appear once per variable it was a candidate for (Variable
    # column) — collapse to one row per gene, taking the first non-NA
    # Acronym/Orthologs across those rows (the two should agree regardless
    # of which variable's row supplied them, since both are gene-level
    # facts, not variable-level ones).
    first_non_na <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA_character_ }
    gs_collapsed <- gs %>%
        group_by(Gene_Name) %>%
        summarise(Acronym = first_non_na(Acronym), Orthologs = first_non_na(Orthologs)) %>%
        ungroup()
    # Acronym preferred; fall back to Orthologs name only when no acronym exists.
    label <- ifelse(!is.na(gs_collapsed$Acronym), gs_collapsed$Acronym,
                    ifelse(!is.na(gs_collapsed$Orthologs), gs_collapsed$Orthologs, NA))
    gene_ortho_label <- setNames(label, gs_collapsed$Gene_Name)
    gene_ortho_label <- gene_ortho_label[!is.na(gene_ortho_label)]
} else {
    print(paste0("WARNING: ", summary_file, " not found — leadingEdge genes will not be annotated with ortholog names"))
}

annotate_leading_edge <- function(genes_str) {
    genes <- strsplit(genes_str, ", ")[[1]]
    labeled <- vapply(genes, function(g) {
        lbl <- if (g %in% names(gene_ortho_label)) gene_ortho_label[[g]] else NA
        if (is.na(lbl)) g else paste0(g, ": ", lbl)
    }, character(1))
    paste(labeled, collapse = ", ")
}
fgsea_res$leadingEdge_orthologs <- vapply(fgsea_res$leadingEdge, annotate_leading_edge, character(1))

outfile <- file.path(output_dir, paste0("GSEA_GO_results_", var_name, tag_suffix, ".tsv"))
write.table(fgsea_res, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)

n_sig <- sum(fgsea_res$padj < 0.05, na.rm = TRUE)
print(paste0("Wrote ", nrow(fgsea_res), " terms with raw p < ", P_MAX,
             " (of ", n_tested, " categories tested) to ", outfile))
print(paste0("Categories tested: ", n_tested,
             " | Min raw p-value: ", format(min(fgsea_res$pval, na.rm = TRUE), digits = 3),
             " | Min FDR (padj): ", format(min(fgsea_res$padj, na.rm = TRUE), digits = 3),
             " | FDR-significant (padj < 0.05): ", n_sig))

# --- Split into per-ontology (BP / MF / CC) result files ---
for (GROUP in c("BP", "MF", "CC")) {
    sub_res <- fgsea_res[!is.na(fgsea_res$ontology) & fgsea_res$ontology == GROUP, ]
    ont_outfile <- file.path(output_dir, paste0("GSEA_GO_results_", var_name, tag_suffix, "_", GROUP, ".tsv"))
    write.table(sub_res, file = ont_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
    print(paste0("Wrote ", nrow(sub_res), " ", GROUP, " terms to ", ont_outfile))
}
n_no_ontology <- sum(is.na(fgsea_res$ontology))
if (n_no_ontology > 0) {
    print(paste0(n_no_ontology, " terms had no ontology (obsolete GO id with no name/ontology match) — ",
                 "included in the combined file but not in any BP/MF/CC split"))
}
