#!/usr/bin/env Rscript --vanilla

# Maps IDs extracted from the SnpEff ANN field (one per line, stdin) to gene IDs
# using the reference GFF, and prints the sorted, deduplicated, dot-suffixed
# gene list to stdout. Needed because some transcript IDs are not numbered
# after their parent gene (e.g. g14399.t1 belongs to gene g4525).
#
# Usage: cat raw_ids.txt | Rscript resolve_gene_ids.R > gene_list.txt

BASE_DIR <- "/private/groups/meyerlab/esteiger"

GFF <- file.path(BASE_DIR, "genome/remove_duplicates_rm_internalstops_longest_isoform.gff")

gff_all <- read.table(GFF, sep = "\t", header = FALSE, quote = "", comment.char = "#",
                      stringsAsFactors = FALSE,
                      col.names = c("scaffold", "source", "feature", "start", "end",
                                    "score", "strand", "frame", "attributes"))

gff_gene_ids <- unique(sub(".*ID=([^;]+).*", "\\1",
                            gff_all$attributes[gff_all$feature == "gene"]))

# tx2gene uses the mRNA ID=/Parent= attributes; the gene_id= attribute can be stale in this GFF.
mrna_rows <- gff_all[gff_all$feature == "mRNA", ]
tx_ids <- sub(".*ID=([^;]+).*", "\\1", mrna_rows$attributes)
gene_ids_col <- sub(".*Parent=([^;]+).*", "\\1", mrna_rows$attributes)
has_both <- grepl("ID=", mrna_rows$attributes) & grepl("Parent=", mrna_rows$attributes)
tx2gene <- setNames(gene_ids_col[has_both], tx_ids[has_both])
tx2gene <- tx2gene[!duplicated(names(tx2gene))]

resolve_one <- function(g) {
    # try the id as-is, then common AUGUSTUS/transdecoder transcript suffixes
    cand <- c(g, paste0(g, ".t1"), paste0(g, ".1"))
    hit <- cand[cand %in% names(tx2gene)]
    if (length(hit) > 0) return(unname(tx2gene[hit[1]]))
    # fall back to a prefix match, for truncated transcript ids like
    # "DEINVI_00006.26586" (real transcript_id "DEINVI_00006.26586.1.p1")
    pref <- grep(paste0("^", gsub("\\.", "\\\\.", g), "\\."), names(tx2gene), value = TRUE)
    if (length(pref) > 0) return(unname(tx2gene[pref[1]]))
    NA_character_
}

raw_ids <- trimws(readLines("stdin"))
raw_ids <- raw_ids[nzchar(raw_ids)]
raw_ids <- sub("\\.$", "", raw_ids)

not_a_gene_id <- !(raw_ids %in% gff_gene_ids)
resolved <- vapply(raw_ids[not_a_gene_id], resolve_one, character(1))
n_resolved <- sum(!is.na(resolved))
raw_ids[not_a_gene_id][!is.na(resolved)] <- resolved[!is.na(resolved)]

if (n_resolved > 0) {
    message(sprintf("resolve_gene_ids.R: resolved %d/%d input ids from transcript to gene id",
                     n_resolved, length(raw_ids)))
}

final <- sort(unique(raw_ids))
cat(paste0(final, "."), sep = "\n")
