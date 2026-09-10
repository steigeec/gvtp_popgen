#!/usr/bin/env Rscript --vanilla

# Cross-variable convergence check for gsea_go_stricter_nominal.R's
# stricter-threshold output (GSEA_GO_stricter_nominal_<var>.tsv). A single
# variable's nominal hit
# is weak evidence on its own (nothing here is FDR-significant, see
# GSEA_GO.R/cluster_gsea_go_results.R); a term or candidate gene that shows
# up as a nominal hit under multiple INDEPENDENT environmental gradients is
# stronger evidence than any one test's p-value, since the three GEMMA runs
# it's drawn from don't share a null distribution.
#
# Reports two levels of convergence:
#   - term-level: the exact same GO term nominal in >=2 variables
#   - gene-level: a leading-edge gene recurring across nominal terms from
#     >=2 different variables (catches convergence GO's granularity misses,
#     e.g. the same gene turning up under two differently-named but related
#     terms in different variables)
#
# Usage: Rscript gsea_go_convergence.R <var1> <var2> [<var3> ...] [--dir <base_dir>]
#   e.g. Rscript gsea_go_convergence.R fog BIO15 BIO18
# Expects each <var>/GSEA_GO_stricter_nominal_<var>.tsv (gsea_go_stricter_nominal.R
# output) under base_dir (default: current directory, one subdirectory per
# variable, matching this pipeline's existing <var>/<var>_file.tsv layout).
# Writes to base_dir:
#   GO_term_convergence_across_variables.tsv
#   GO_gene_convergence_across_variables.tsv

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
dir_idx <- which(args == "--dir")
base_dir <- if (length(dir_idx) == 1 && length(args) >= dir_idx + 1) args[dir_idx + 1] else "."
var_names <- if (length(dir_idx) == 1) args[-c(dir_idx, dir_idx + 1)] else args
if (length(var_names) < 2) {
    stop("Usage: Rscript gsea_go_convergence.R <var1> <var2> [<var3> ...] [--dir <base_dir>]")
}

# --- Load each variable's stricter-threshold nominal table ---
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

# --- Term-level convergence: same GO id (pathway) nominal in >=2 variables ---
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
print(paste0(nrow(term_conv), " GO terms nominal in >=2 of ", length(var_names),
             " variables -- wrote ", term_out))

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
