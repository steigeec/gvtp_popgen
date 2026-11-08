#!/usr/bin/env Rscript --vanilla

# Dump AGI gene_id -> GO term pairs from org.At.tair.db, for
# build_tair_ortholog_go.py. Uses direct annotations (GO, not GOALL) and
# drops IEA and ND evidence (ND rows point only at root terms).
# Usage: Rscript dump_tair_go.R [output.tsv]

BASE_DIR <- "/private/groups/meyerlab/esteiger"

suppressMessages({
    library(org.At.tair.db)
    library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1) args[1] else
    file.path(BASE_DIR, "genome/annotation/tair_geneid_go.tsv")

all_tair <- keys(org.At.tair.db, keytype = "TAIR")
df <- suppressWarnings(AnnotationDbi::select(
    org.At.tair.db, keys = all_tair, keytype = "TAIR",
    columns = c("GO", "EVIDENCE")
))

n_raw <- nrow(df)
df <- df[!is.na(df$GO) & !(df$EVIDENCE %in% c("IEA", "ND")), c("TAIR", "GO")]
df <- unique(df)

write.table(df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = c("gene_id", "go_term"))

print(paste0("org.At.tair.db: ", n_raw, " raw TAIR/GO/evidence rows -> ",
             nrow(df), " unique gene-GO pairs (IEA + ND evidence excluded) ",
             "across ", length(unique(df$TAIR)), " AT genes"))
print(paste0("Wrote ", out_file))
