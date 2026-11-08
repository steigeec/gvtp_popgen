#!/usr/bin/env Rscript --vanilla

# Selects GSEA_GO.R terms at a stricter nominal threshold (raw p < p_max) than
# cluster_gsea_go_results.R, then applies the same rrvgo redundancy reduction
# (per ontology, Wang similarity, threshold 0.7). Output is plotted by
# gsea_go_lollipop_3panel.R.
#
# Usage: Rscript gsea_go_stricter_nominal.R <var_name> [output_dir] [p_max] [tag]
#   e.g. Rscript gsea_go_stricter_nominal.R fog fog 0.005
# Expects GSEA_GO_results_<var_name>_{BP,MF,CC}.tsv (GSEA_GO.R output) in
# output_dir (default: <var_name>).
# Writes:
#   <output_dir>/GSEA_GO_stricter_nominal_<var_name>.tsv

suppressMessages({
    library(dplyr)
    library(rrvgo)
})

SIM_THRESHOLD <- 0.7  # same as cluster_gsea_go_results.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript gsea_go_stricter_nominal.R <var_name> [output_dir] [p_max] [tag]")
}
var_name   <- args[1]
output_dir <- if (length(args) >= 2) args[2] else var_name
P_MAX      <- if (length(args) >= 3) as.numeric(args[3]) else 0.005
# tag must match the output_tag used for GSEA_GO.R (e.g. "sunflower").
tag        <- if (length(args) >= 4) args[4] else ""
tag_suffix <- if (tag != "") paste0("_", tag) else ""

# --- Load and combine all three ontologies ---
all_res <- list()
for (GROUP in c("BP", "MF", "CC")) {
    infile <- file.path(output_dir, paste0("GSEA_GO_results_", var_name, tag_suffix, "_", GROUP, ".tsv"))
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

# --- Nominal selection, then rrvgo reduction per ontology ---
reduced_list <- list()
n_pooled <- 0
for (GROUP in names(all_res)) {
    sub_res <- combined[combined$ontology == GROUP, ]
    sub_res <- sub_res[sub_res$pval < P_MAX, ]
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
print(paste0("Nominal hits (raw p < ", P_MAX, "), pooled across ontologies: ",
             n_pooled, " -> ", nrow(reduced), " after rrvgo redundancy reduction"))

out_cols <- c("pathway", "term", "ontology", "pval", "padj", "NES", "size",
              "leadingEdge", "leadingEdge_orthologs")
out_cols <- intersect(out_cols, names(reduced))
reduced_out <- reduced[, out_cols]
outfile <- file.path(output_dir, paste0("GSEA_GO_stricter_nominal_", var_name, tag_suffix, ".tsv"))
write.table(reduced_out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
print(paste0("Wrote ", nrow(reduced_out), " nominal hits (raw p < ", P_MAX,
             ", NOT FDR-significant) to ", outfile))
