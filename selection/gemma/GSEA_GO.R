#!/usr/bin/env Rscript --vanilla

# Preranked GO enrichment (fgsea) with all genes ranked by their per-gene
# ACAT-combined GEMMA p-value (build_gene_gwas_scores.py).
#
# Usage: Rscript GSEA_GO.R <var_name> [output_dir] [go_annot_file] [output_tag]
#   e.g. Rscript GSEA_GO.R fog fog
# Expects <var_name>_gene_gwas_scores.tsv (build_gene_gwas_scores.py output)
# in output_dir (default: <var_name>).

BASE_DIR <- "/private/groups/meyerlab/esteiger"

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
# go_annot_file: alternate gene-GO annotation source (default: combined annotation).
# output_tag: appended to output filenames, e.g. GSEA_GO_results_<var>_<tag>.tsv
go_annot_file <- if (length(args) >= 3) args[3] else
    file.path(BASE_DIR, "genome/annotated_sequences_combined_geneid_go.tsv")
output_tag  <- if (length(args) >= 4) args[4] else ""
tag_suffix  <- if (output_tag != "") paste0("_", output_tag) else ""

scores_file <- file.path(output_dir, paste0(var_name, "_gene_gwas_scores.tsv"))

# --- Load GO annotation ---
go_annot_raw <- read.table(go_annot_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
names(go_annot_raw) <- c("gene", "go")
go_annot_raw$gene <- trimws(go_annot_raw$gene)
go_annot_raw <- unique(go_annot_raw[, c("gene", "go")])
print(paste0("Loaded ", nrow(go_annot_raw), " gene-GO pairs across ",
             length(unique(go_annot_raw$gene)), " genes"))

# Remap obsolete GO ids to their current successor (unambiguous replaced_by chains;
# table from build_go_ancestor_map.py). Obsolete ids cannot be propagated to ancestors.
remap_file <- file.path(BASE_DIR, "genome/go_obsolete_remap.tsv")
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

# Drop terms implausible in plants before propagation and testing: union of
# build_nonplant_go_exclusion_list.py (anchor-term heuristic) and
# build_taxon_constraint_go_exclusion_list.py (GO never_in_taxon Viridiplantae).
excl_file       <- file.path(BASE_DIR, "genome/go_terms_excluded_nonplant.txt")
excl_file_taxon <- file.path(BASE_DIR, "genome/go_terms_excluded_taxon_constraint.txt")
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

# Also drop generic housekeeping terms (translation, proteasome, core metabolism,
# ontology roots, etc.; build_housekeeping_go_exclusion_list.py) before testing.
housekeeping_file <- file.path(BASE_DIR, "genome/go_terms_excluded_housekeeping.txt")
if (file.exists(housekeeping_file)) {
    housekeeping_ids <- read.table(housekeeping_file, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")[[1]]
    n_before <- nrow(go_annot_raw)
    go_annot_raw <- go_annot_raw[!(go_annot_raw$go %in% housekeeping_ids), ]
    print(paste0("Housekeeping term filter: ", n_before, " -> ", nrow(go_annot_raw), " direct gene-GO pairs"))
} else {
    print(paste0("WARNING: housekeeping exclusion list not found at ", housekeeping_file, " — skipping that filter"))
}

# --- Propagate annotations up the GO hierarchy (true path rule) ---
# Ancestor map from a current go-basic.obo (build_go_ancestor_map.py), not GO.db
anc_tab <- read.table(file.path(BASE_DIR, "genome/go_ancestor_map.tsv"),
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

# Re-apply the housekeeping filter: propagation adds back root and near-root ancestors
if (file.exists(housekeeping_file)) {
    n_before <- nrow(go_annot_prop)
    go_annot_prop <- go_annot_prop[!(go_annot_prop$go %in% housekeeping_ids), ]
    print(paste0("Housekeeping term filter (post-propagation): ", n_before, " -> ", nrow(go_annot_prop), " gene-term pairs"))
}
# Same re-pass for non-plant terms
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
# Sized against the full annotated gene set
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

# Rank by snp_count_adjusted_score (build_gene_gwas_scores.py): -log10(ACAT p) with a
# LOWESS fit on log SNP count removed. The residual can be negative, so it is
# converted to ranks (1..N) for scoreType = "pos".
stats <- rank(scores$snp_count_adjusted_score)
names(stats) <- scores$gene_id
stats <- sort(stats, decreasing = TRUE)

# --- Run fgsea ---
# scoreType = "pos": the ranking statistic is one-tailed (no "down" direction)
set.seed(1)
fgsea_res <- fgsea(pathways = pathways, stats = stats, minSize = MIN_CAT_SIZE,
                    scoreType = "pos", eps = 0)
fgsea_res <- as.data.frame(fgsea_res)  # data.table column indexing differs from data.frame
fgsea_res <- fgsea_res[order(fgsea_res$pval), ]
fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, paste, collapse = ", ")
n_tested <- nrow(fgsea_res)

# --- Human-readable term names + ontology (BP/MF/CC) ---
# From go_term_names.tsv (build_go_ancestor_map.py), which is more current than GO.db;
# obsolete GO ids get NA.
term_name_tab <- read.table(file.path(BASE_DIR, "genome/go_term_names.tsv"),
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "",
                            comment.char = "")  # some term names contain "#" (e.g. "ascr#2 binding")
match_idx <- match(fgsea_res$pathway, term_name_tab$id)
fgsea_res$term <- term_name_tab$name[match_idx]
fgsea_res$ontology <- term_name_tab$ontology[match_idx]
fgsea_res <- fgsea_res[, c("pathway", "term", "ontology", setdiff(names(fgsea_res), c("pathway", "term", "ontology")))]

# --- Report only nominal hits (raw p < P_MAX) ---
fgsea_res <- fgsea_res[!is.na(fgsea_res$pval) & fgsea_res$pval < P_MAX, ]

# --- Label leadingEdge genes with acronym or ortholog name from gene_summary_all.tsv ---
# That file covers only candidate genes, so most leadingEdge genes keep a bare gene id.
summary_file <- file.path(BASE_DIR, "selection/gemma/gene_summary_all.tsv")
gene_ortho_label <- character(0)
if (file.exists(summary_file)) {
    gs <- read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                     quote = "", comment.char = "")
    # One row per gene: first non-NA Acronym/Orthologs across variables
    first_non_na <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA_character_ }
    gs_collapsed <- gs %>%
        group_by(Gene_Name) %>%
        summarise(Acronym = first_non_na(Acronym), Orthologs = first_non_na(Orthologs)) %>%
        ungroup()
    # Acronym preferred, else Orthologs name
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
