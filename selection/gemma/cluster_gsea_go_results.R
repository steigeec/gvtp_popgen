#!/usr/bin/env Rscript --vanilla

# Reduces GSEA_GO.R's raw per-term output (GSEA_GO_results_<var>_{BP,MF,CC}.tsv
# -- hundreds of GO terms, many nested parent/child chains driven by
# substantially the same underlying gene(s) via propagation) down to one
# representative term per semantically-distinct cluster, using the same
# rrvgo approach top_nominal_GO_hits.R already applies to the threshold-
# based (goseq) track. GSEA_GO_results is already clean of non-plant and
# housekeeping terms (GSEA_GO.R drops both BEFORE testing, not just from
# the report), so this script only needs to select + cluster, not re-filter.
#
# No FDR-significant terms exist in this dataset (min padj ~0.37-0.67
# across all three ontologies and all three environmental variables, per
# GSEA_GO.R's own printed summary) -- this is a nominal/hypothesis-
# generating selection (raw p < P_NOMINAL_MAX), explicitly not a
# substitute for the FDR-based result, matching the same caveat
# top_nominal_GO_hits.R carries for its own track.
#
# Combines all three ontologies into ONE output table (not three separate
# ones, and not BP-only) with an Ontology column -- restricting to BP would
# have hidden real hits that only show up in MF/CC (e.g. this dataset's
# "photoreceptor activity" MF term, which is one of the most mechanistically
# specific hits in the whole set, alongside the BP-level light-signaling
# terms it corroborates).
#
# Usage: Rscript cluster_gsea_go_results.R <var_name> [output_dir] [p_max]
#   e.g. Rscript cluster_gsea_go_results.R fog fog
# Expects GSEA_GO_results_<var_name>_{BP,MF,CC}.tsv (GSEA_GO.R output) in
# output_dir (default: same directory this is run from, matching
# GSEA_GO.R's own default output location).

suppressMessages({
    library(dplyr)
    library(rrvgo)
})

SIM_THRESHOLD <- 0.7  # same clustering threshold top_nominal_GO_hits.R uses

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript cluster_gsea_go_results.R <var_name> [output_dir] [p_max]")
}
var_name      <- args[1]
output_dir    <- if (length(args) >= 2) args[2] else var_name
P_NOMINAL_MAX <- if (length(args) >= 3) as.numeric(args[3]) else 0.05

# --- Load + combine all three ontologies ---
all_res <- list()
for (GROUP in c("BP", "MF", "CC")) {
    infile <- file.path(output_dir, paste0("GSEA_GO_results_", var_name, "_", GROUP, ".tsv"))
    if (!file.exists(infile)) {
        print(paste0("WARNING: ", infile, " not found -- skipping ", GROUP))
        next
    }
    df <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
    df$ontology <- GROUP
    all_res[[GROUP]] <- df
}
if (length(all_res) == 0) {
    stop("No GSEA_GO_results_*.tsv files found for ", var_name, " in ", output_dir)
}
combined <- do.call(rbind, all_res)
combined <- combined[!is.na(combined$pval), ]
print(paste0(nrow(combined), " terms loaded across ", length(all_res), " ontologies"))

# --- Nominal selection, per ontology (rrvgo's calculateSimMatrix requires a
# single ont at a time -- can't mix BP/MF/CC terms in one similarity matrix) ---
reduced_list <- list()
n_pooled <- 0
for (GROUP in names(all_res)) {
    sub_res <- combined[combined$ontology == GROUP, ]
    sub_res <- sub_res[sub_res$pval < P_NOMINAL_MAX, ]
    n_pooled <- n_pooled + nrow(sub_res)
    if (nrow(sub_res) >= 3) {
        simMatrix <- calculateSimMatrix(sub_res$pathway, orgdb = "org.At.tair.db",
                                        ont = GROUP, method = "Wang")
        reducedTerms <- reduceSimMatrix(simMatrix,
                                        scores = setNames(-log10(sub_res$pval), sub_res$pathway),
                                        threshold = SIM_THRESHOLD, orgdb = "org.At.tair.db")
        reduced_go <- reducedTerms[reducedTerms$parent == reducedTerms$go, ]
        sub_res <- sub_res[sub_res$pathway %in% reduced_go$go, ]
    }
    reduced_list[[GROUP]] <- sub_res
}
reduced <- do.call(rbind, reduced_list)
reduced <- reduced[order(reduced$pval), ]
print(paste0("Nominal hits (raw p < ", P_NOMINAL_MAX, "), pooled across ontologies: ",
             n_pooled, " -> ", nrow(reduced), " after rrvgo redundancy reduction"))

# --- Write combined, ranked table ---
out_cols <- c("pathway", "term", "ontology", "pval", "padj", "NES", "size",
              "leadingEdge", "leadingEdge_orthologs")
out_cols <- intersect(out_cols, names(reduced))
reduced <- reduced[, out_cols]

outfile <- file.path(output_dir, paste0("GSEA_GO_top_nominal_", var_name, ".tsv"))
write.table(reduced, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
print(paste0("Wrote ", nrow(reduced), " nominal hits (raw p < ", P_NOMINAL_MAX,
             ", NOT FDR-significant -- see module docstring) to ", outfile))
print(paste0("By ontology: ", paste(names(table(reduced$ontology)), table(reduced$ontology),
                                    sep = "=", collapse = ", ")))
