#!/usr/bin/env Rscript --vanilla

# Cross-variable convergence of nominal GSEA GO results
# (GSEA_GO_stricter_nominal_<var>.tsv from gsea_go_stricter_nominal.R).
# Reports:
#   - term-level, exact match: the same GO term nominal in >=2 variables;
#   - term-level near-miss candidates: pairs of different GO terms, one from
#     each of two variables, with Wang similarity >= NEAR_MISS_THRESHOLD
#     (same ontology). Separate rrvgo reductions can keep different sibling
#     terms for the same cluster. These are for manual review, not merged;
#   - gene-level: a leading-edge gene recurring across nominal terms from
#     >=2 variables.
#
# Usage: Rscript gsea_go_convergence.R <var1> <var2> [<var3> ...] [--dir <base_dir>]
#   e.g. Rscript gsea_go_convergence.R fog BIO15 BIO18
# Expects each <var>/GSEA_GO_stricter_nominal_<var>.tsv (gsea_go_stricter_nominal.R
# output) under base_dir (default: current directory, one subdirectory per variable).
# Writes to base_dir:
#   GO_term_convergence_across_variables.tsv   (exact-match, primary)
#   GO_term_near_miss_candidates.tsv           (semantic near-misses, for manual review)
#   GO_gene_convergence_across_variables.tsv

suppressMessages({
    library(dplyr)
    library(GOSemSim)
})

NEAR_MISS_THRESHOLD <- 0.8  # Wang similarity

args <- commandArgs(trailingOnly = TRUE)
dir_idx <- which(args == "--dir")
base_dir <- if (length(dir_idx) == 1 && length(args) >= dir_idx + 1) args[dir_idx + 1] else "."
var_names <- if (length(dir_idx) == 1) args[-c(dir_idx, dir_idx + 1)] else args
if (length(var_names) < 2) {
    stop("Usage: Rscript gsea_go_convergence.R <var1> <var2> [<var3> ...] [--dir <base_dir>]")
}

# --- Load each variable's rrvgo-reduced nominal table ---
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

# --- Term-level convergence (exact match): same GO id nominal in >=2 variables ---
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

# --- Term-level near-miss candidates: different GO ids, same ontology,
# high pairwise semantic similarity ---
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
# terms from >=2 variables ---
# Terms larger than MAX_TERM_SIZE (e.g. "plasma membrane") are excluded from this
# step only, since very large terms share genes with almost anything.
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
