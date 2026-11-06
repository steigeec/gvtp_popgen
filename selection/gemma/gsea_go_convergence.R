#!/usr/bin/env Rscript --vanilla

# Cross-variable convergence check for gsea_go_stricter_nominal.R's
# stricter-threshold output (GSEA_GO_stricter_nominal_<var>.tsv). A single
# variable's nominal hit is weak evidence on its own (nothing here is
# FDR-significant, see GSEA_GO.R/cluster_gsea_go_results.R); a term or
# candidate gene that shows up as a nominal hit under multiple INDEPENDENT
# environmental gradients is stronger evidence than any one test's p-value,
# since the three GEMMA runs it's drawn from don't share a null
# distribution.
#
# Reports THREE things:
#   - term-level (exact match): the exact same GO term nominal in >=2
#     variables' own already-reduced lists. This is the primary,
#     highest-confidence convergence result.
#   - term-level (near-miss candidates): pairs of DIFFERENT GO ids, one
#     from each of two variables' final lists, whose GO semantic
#     similarity (Wang, same ontology) is >= NEAR_MISS_THRESHOLD. This
#     exists because gsea_go_stricter_nominal.R's rrvgo reduction runs
#     independently per variable, and its greedy highest-scoring-member
#     rule can pick a DIFFERENT sibling term as the "representative" in
#     each variable even when the same semantic cluster is enriched in
#     both -- confirmed on this dataset: wax biosynthesis survives as
#     GO:1904278 "positive regulation of wax biosynthetic process" in
#     fog's reduction but GO:1904276 "regulation of wax biosynthetic
#     process" in BIO18's, so the exact-match check above alone would
#     silently miss it. This near-miss list is intentionally NOT
#     auto-merged into the exact-match table: pairwise GO similarity
#     doesn't cleanly separate "same process" from "related but
#     meaningfully different" (empirically, a probably-genuine pair like
#     diterpene phytoalexin biosynthetic/metabolic process scored LOWER,
#     0.421, than a probably-different pair like protein arginylation vs.
#     glutathionylation, 0.621) -- so these candidates need a human
#     decision, not an automatic one. Also deliberately pairwise rather
#     than clustering-based: clustering with transitive chaining is what
#     caused the earlier over-merging (e.g. "protein arginylation" and
#     "positive regulation of mRNA catabolic process" ending up in one
#     cluster at rrvgo threshold 0.7 despite a direct similarity of only
#     0.214) -- a direct pairwise check has no transitive chain to ride.
#   - gene-level: a leading-edge gene recurring across nominal terms from
#     >=2 different variables (catches convergence GO's granularity misses
#     entirely, e.g. the same gene turning up under two differently-named
#     but related terms in different variables)
#
# Usage: Rscript gsea_go_convergence.R <var1> <var2> [<var3> ...] [--dir <base_dir>]
#   e.g. Rscript gsea_go_convergence.R fog BIO15 BIO18
# Expects each <var>/GSEA_GO_stricter_nominal_<var>.tsv (gsea_go_stricter_nominal.R
# output) under base_dir (default: current directory, one subdirectory per
# variable, matching this pipeline's existing <var>/<var>_file.tsv layout).
# Writes to base_dir:
#   GO_term_convergence_across_variables.tsv   (exact-match, primary)
#   GO_term_near_miss_candidates.tsv           (semantic near-misses, for manual review)
#   GO_gene_convergence_across_variables.tsv

suppressMessages({
    library(dplyr)
    library(GOSemSim)
})

NEAR_MISS_THRESHOLD <- 0.8  # Wang similarity; see header note on why this is
                            # a manual-review candidate list, not an
                            # auto-merge -- 0.8 catches the confirmed wax
                            # biosynthesis case (0.827) while excluding a
                            # probably-different pair scoring 0.621.

args <- commandArgs(trailingOnly = TRUE)
dir_idx <- which(args == "--dir")
base_dir <- if (length(dir_idx) == 1 && length(args) >= dir_idx + 1) args[dir_idx + 1] else "."
var_names <- if (length(dir_idx) == 1) args[-c(dir_idx, dir_idx + 1)] else args
if (length(var_names) < 2) {
    stop("Usage: Rscript gsea_go_convergence.R <var1> <var2> [<var3> ...] [--dir <base_dir>]")
}

# --- Load each variable's stricter-threshold nominal table (each
# variable's own, independently rrvgo-reduced top-term list -- the exact
# same file gsea_go_lollipop_3panel.R plots) ---
all_terms <- list()
for (v in var_names) {
    infile <- file.path(base_dir, v, paste0("GSEA_GO_stricter_nominal_", v, ".tsv"))
    if (!file.exists(infile)) {
        stop("Missing ", infile, " -- run gsea_go_stricter_nominal.R for '", v, "' first.")
    }
    df <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
    df$variable <- v
    all_terms[[v]] <- df
}
combined <- bind_rows(all_terms)
print(paste0("Loaded ", nrow(combined), " nominal term-hits across ", length(var_names), " variables"))

# --- Term-level convergence (primary, exact match): same GO id (pathway)
# nominal in >=2 variables ---
term_conv <- combined %>%
    group_by(pathway, term, ontology) %>%
    summarise(
        n_variables = n_distinct(variable),
        variables   = paste(sort(unique(variable)), collapse = ", "),
        min_pval    = min(pval),
        pvals_by_variable = paste(sprintf("%s=%.4g", variable, pval), collapse = "; "),
        .groups = "drop"
    ) %>%
    filter(n_variables >= 2) %>%
    arrange(desc(n_variables), min_pval)

term_out <- file.path(base_dir, "GO_term_convergence_across_variables.tsv")
write.table(term_conv, file = term_out, sep = "\t", quote = FALSE, row.names = FALSE)
print(paste0(nrow(term_conv), " GO terms (exact match) nominal in >=2 of ", length(var_names),
             " variables -- wrote ", term_out))

# --- Term-level near-miss candidates: different GO id, high pairwise
# semantic similarity, one term from each of two variables. Computed only
# between variable-pairs' terms that AREN'T already an exact match above,
# and only within the same ontology (cross-ontology similarity isn't
# meaningful for GOSemSim's annotation-based methods). ---
near_miss_rows <- list()
ontologies_present <- intersect(c("BP", "MF", "CC"), unique(combined$ontology))
sem_data <- setNames(lapply(ontologies_present, function(o) godata("org.At.tair.db", ont = o)),
                      ontologies_present)

var_pairs <- combn(var_names, 2, simplify = FALSE)
for (pair in var_pairs) {
    v1 <- combined[combined$variable == pair[1], ]
    v2 <- combined[combined$variable == pair[2], ]
    for (ont in ontologies_present) {
        t1 <- v1[v1$ontology == ont, ]
        t2 <- v2[v2$ontology == ont, ]
        if (nrow(t1) == 0 || nrow(t2) == 0) next
        for (i in seq_len(nrow(t1))) {
            for (j in seq_len(nrow(t2))) {
                if (t1$pathway[i] == t2$pathway[j]) next  # already in the exact-match table
                sim <- tryCatch(
                    goSim(t1$pathway[i], t2$pathway[j], semData = sem_data[[ont]], measure = "Wang"),
                    error = function(e) NA_real_
                )
                if (!is.na(sim) && sim >= NEAR_MISS_THRESHOLD) {
                    near_miss_rows[[length(near_miss_rows) + 1]] <- data.frame(
                        ontology = ont, similarity = round(sim, 3),
                        variable_1 = pair[1], pathway_1 = t1$pathway[i], term_1 = t1$term[i], pval_1 = t1$pval[i],
                        variable_2 = pair[2], pathway_2 = t2$pathway[j], term_2 = t2$term[j], pval_2 = t2$pval[j],
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }
}
near_miss <- if (length(near_miss_rows) > 0) {
    bind_rows(near_miss_rows) %>% arrange(desc(similarity))
} else {
    data.frame(ontology = character(), similarity = numeric(),
               variable_1 = character(), pathway_1 = character(), term_1 = character(), pval_1 = numeric(),
               variable_2 = character(), pathway_2 = character(), term_2 = character(), pval_2 = numeric())
}

near_miss_out <- file.path(base_dir, "GO_term_near_miss_candidates.tsv")
write.table(near_miss, file = near_miss_out, sep = "\t", quote = FALSE, row.names = FALSE)
print(paste0(nrow(near_miss), " near-miss candidate pairs (similarity >= ", NEAR_MISS_THRESHOLD,
             ", different GO ids, NOT auto-merged -- manual review) -- wrote ", near_miss_out))

# --- Gene-level convergence: a leading-edge gene recurring across nominal
# terms from >=2 different variables (regardless of which term) ---
# Terms above MAX_TERM_SIZE are excluded from this step only (not from the
# term-level convergence above): a near-generic category like "plasma
# membrane" (size ~4091, a large fraction of the annotated genome here) will
# trivially "converge" with almost any gene in almost any variable just by
# being huge, not because of any real shared signal -- confirmed empirically
# on this dataset, where excluding it alone dropped 141 of 249 gene-level
# convergence rows (57%). fgsea itself already corrects for set size when
# testing a term's own enrichment; this cap only protects the separate,
# unweighted gene-overlap-counting logic below, which has no such correction.
MAX_TERM_SIZE <- 200

n_before_size_filter <- nrow(combined)
combined_for_genes <- combined[combined$size <= MAX_TERM_SIZE, ]
print(paste0("Gene-level convergence: excluding terms with size > ", MAX_TERM_SIZE,
             " (", n_before_size_filter - nrow(combined_for_genes), " of ", n_before_size_filter,
             " term-hits dropped, e.g. generic categories like 'plasma membrane')"))

gene_rows <- list()
for (i in seq_len(nrow(combined_for_genes))) {
    genes <- unique(trimws(strsplit(combined_for_genes$leadingEdge[i], ",\\s*")[[1]]))
    if (length(genes) == 0) next
    gene_rows[[length(gene_rows) + 1]] <- data.frame(
        gene     = genes,
        variable = combined_for_genes$variable[i],
        term     = combined_for_genes$term[i],
        ontology = combined_for_genes$ontology[i],
        pval     = combined_for_genes$pval[i],
        stringsAsFactors = FALSE
    )
}
gene_long <- bind_rows(gene_rows)

gene_conv <- gene_long %>%
    group_by(gene) %>%
    summarise(
        n_variables = n_distinct(variable),
        variables   = paste(sort(unique(variable)), collapse = ", "),
        n_terms     = n_distinct(term),
        terms       = paste(sort(unique(term)), collapse = "; "),
        .groups = "drop"
    ) %>%
    filter(n_variables >= 2) %>%
    arrange(desc(n_variables), desc(n_terms))

gene_out <- file.path(base_dir, "GO_gene_convergence_across_variables.tsv")
write.table(gene_conv, file = gene_out, sep = "\t", quote = FALSE, row.names = FALSE)
print(paste0(nrow(gene_conv), " leading-edge genes recurring across nominal terms from >=2 of ",
             length(var_names), " variables -- wrote ", gene_out))
