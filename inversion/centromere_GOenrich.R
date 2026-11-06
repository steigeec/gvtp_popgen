#!/usr/bin/env Rscript --vanilla

# GO enrichment of the genes in the chromosome 6 centromeric low-recombination
# block, which includes the core fog-outlier block (49.0-49.2 Mb; the
# tightly linked haplotype carrying GNS1 and SNF1) and the centromeric
# low-recombination block (44.6-48.9 Mb; moderate LD with the fog haplotype,
# containing the LINE-rich centromere at ~45.0-48.5 Mb). See
# inversions/results.txt for how these limits were defined.
#
# Binary test (gene in the block vs not), NOT a GSEA: target = every gene
# (GFF "gene" feature) on BLOCK_CHROM overlapping a BLOCK interval; universe
# = all GO-annotated D. increscens subsp. villosa genes with a GFF length.
# Gene-length-corrected with goseq (Wallenius approximation): a probability
# weighting function on gene length (genomic span) is fitted before testing,
# so the result is not driven by the block having unusually long or short
# genes.
#
# GO annotation and curation follow gvtp_popgen/selection/gemma/GSEA_GO.R
# exactly (the curated annotation used for the GWAS GSEA):
#   - combined gene-level GO annotation (eggNOG-mapper + harvest_blast_go.py
#     NCBI recovery; genome/annotated_sequences_combined_geneid_go.tsv)
#   - obsolete GO ids remapped to their current successor (go_obsolete_remap.tsv)
#   - non-plant-plausible terms removed: union of the heuristic list and the
#     GO Consortium never_in_taxon-Viridiplantae list (before and after
#     propagation)
#   - generic housekeeping terms removed (before and after propagation)
#   - true-path propagation from a fresh go-basic.obo (go_ancestor_map.tsv)
#   - term names and ontology from go_term_names.tsv (same obo), not GO.db
#   - categories with < MIN_CAT_SIZE universe genes dropped
# One deliberate difference from GSEA_GO.R: obsolete ids with no successor
# (no name/ontology in the fresh obo) are dropped BEFORE testing, so they
# don't spend the FDR budget (GSEA_GO.R tests them and leaves them out of
# its per-ontology files afterwards).
# goseq's gene-length correction is from GOenrich_goseq.R. BH FDR across BP,
# MF and CC together; rrvgo (Wang similarity) reduces redundancy among
# FDR-significant terms.
#
# Note: the block is gene-poor, so the target set is small and power is
# limited; the full results table (all tested terms) is written too, so
# nominal patterns can be inspected even if nothing passes FDR.
#
# Usage: Rscript centromere_GOenrich.R [output_dir]

.libPaths(c("/private/groups/meyerlab/esteiger/software/R_libs", .libPaths()))
if (getRversion() < "4.4.0") `%||%` <- function(x, y) if (is.null(x)) y else x

suppressMessages({
    library(goseq)
    library(dplyr)
    library(rrvgo)
    library(GO.db)
})

BLOCK_CHROM <- "DEINVI_00006"
BLOCK <- data.frame(name  = c("centromeric low-recombination block", "core fog-outlier block"),
                    start = c(44.6e6, 49.0e6),
                    end   = c(48.9e6, 49.2e6))
MIN_CAT_SIZE <- 3

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else "/private/groups/meyerlab/esteiger/inversions/centromere/GO"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
file_label <- "chr6_centromeric_block"

GFF       <- "/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"
GO_ANNOT  <- "/private/groups/meyerlab/esteiger/genome/annotated_sequences_combined_geneid_go.tsv"
REMAP     <- "/private/groups/meyerlab/esteiger/genome/go_obsolete_remap.tsv"
EXCL      <- c("/private/groups/meyerlab/esteiger/genome/go_terms_excluded_nonplant.txt",
               "/private/groups/meyerlab/esteiger/genome/go_terms_excluded_taxon_constraint.txt")
HOUSEKEEPING <- "/private/groups/meyerlab/esteiger/genome/go_terms_excluded_housekeeping.txt"
TERM_NAMES   <- "/private/groups/meyerlab/esteiger/genome/go_term_names.tsv"
ANCESTORS <- "/private/groups/meyerlab/esteiger/genome/go_ancestor_map.tsv"


# --- GO annotation and curation (as in GSEA_GO.R) ---
go_annot_raw <- read.table(GO_ANNOT, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
names(go_annot_raw) <- c("gene", "go")
go_annot_raw$gene <- trimws(go_annot_raw$gene)
go_annot_raw <- unique(go_annot_raw[, c("gene", "go")])
print(paste0("Loaded ", nrow(go_annot_raw), " gene-GO pairs across ", length(unique(go_annot_raw$gene)), " genes"))

remap_tab <- read.table(REMAP, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
remap_vec <- setNames(remap_tab$new_id, remap_tab$old_id)
go_annot_raw$go <- ifelse(go_annot_raw$go %in% names(remap_vec), remap_vec[go_annot_raw$go], go_annot_raw$go)
go_annot_raw <- unique(go_annot_raw)

read_ids <- function(f) read.table(f, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")[[1]]
excl_ids <- unique(unlist(lapply(EXCL, read_ids)))
housekeeping_ids <- read_ids(HOUSEKEEPING)
n0 <- nrow(go_annot_raw)
go_annot_raw <- go_annot_raw[!(go_annot_raw$go %in% c(excl_ids, housekeeping_ids)), ]
print(paste0("Non-plant (", length(excl_ids), " ids) + housekeeping (", length(housekeeping_ids),
             " ids) filters: ", n0, " -> ", nrow(go_annot_raw), " direct gene-GO pairs"))

anc_tab <- read.table(ANCESTORS, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
anc_map <- split(anc_tab$ancestor, anc_tab$term)
unique_terms <- unique(go_annot_raw$go)
anc_list <- lapply(unique_terms, function(t) {
    a <- anc_map[[t]]
    if (is.null(a)) return(t)
    unique(c(t, a[grepl("^GO:", a)]))
})
names(anc_list) <- unique_terms
expanded <- anc_list[go_annot_raw$go]
go_annot_prop <- unique(data.frame(gene = rep(go_annot_raw$gene, lengths(expanded)),
                                   go = unlist(expanded, use.names = FALSE), stringsAsFactors = FALSE))
print(paste0("GO propagation: ", nrow(go_annot_raw), " direct gene-term pairs -> ", nrow(go_annot_prop)))
# propagation can re-add housekeeping / non-plant terms as ancestors: second pass
go_annot_prop <- go_annot_prop[!(go_annot_prop$go %in% c(excl_ids, housekeeping_ids)), ]
# names / ontology from the fresh obo; obsolete ids with no successor dropped
term_tab <- read.table(TERM_NAMES, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "",
                       comment.char = "")   # some term names contain "#" (as in GSEA_GO.R)
n_obs <- length(setdiff(unique(go_annot_prop$go), term_tab$id))
go_annot_prop <- go_annot_prop[go_annot_prop$go %in% term_tab$id, ]
print(paste0("After post-propagation filters and dropping ", n_obs, " obsolete ids with no successor: ",
             nrow(go_annot_prop), " gene-term pairs"))
gene2GO <- split(go_annot_prop$go, go_annot_prop$gene)
all_genes <- names(gene2GO)


# --- Genes, lengths and positions from the GFF ---
gff <- read.table(GFF, sep = "\t", header = FALSE, quote = "", comment.char = "#", stringsAsFactors = FALSE,
                  col.names = c("scaffold", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"))
genes <- gff[gff$feature == "gene", ]
genes$gene_id <- sub(".*gene_id=([^;]+).*", "\\1", genes$attributes)
gene_length <- setNames(genes$end - genes$start + 1, genes$gene_id)


# --- Target set: genes overlapping the block ---
in_block <- rep(FALSE, nrow(genes)); block_name <- rep(NA_character_, nrow(genes))
for (i in seq_len(nrow(BLOCK))) {
    hit <- genes$scaffold == BLOCK_CHROM & genes$end >= BLOCK$start[i] & genes$start <= BLOCK$end[i]
    in_block <- in_block | hit
    block_name[hit] <- BLOCK$name[i]
}
block_genes <- genes[in_block, c("gene_id", "scaffold", "start", "end")]
block_genes$part <- block_name[in_block]
block_genes$length <- block_genes$end - block_genes$start + 1
block_genes$has_GO <- block_genes$gene_id %in% all_genes
block_genes$n_GO_direct <- sapply(block_genes$gene_id, function(g) sum(go_annot_raw$gene == g))

# Ortholog names (descriptive only; not used in the test):
#  - best_blast_hit / best_blast_description: top hit (by bitscore) of the gene's
#    annotated transcript in the genome-wide blastp (BLASTP_OUT)
#  - arabidopsis_best_hit / arabidopsis_same_orthogroup: best Arabidopsis hit in
#    GENESPACE's all-vs-all blast, and whether OrthoFinder put the two in one
#    orthogroup (a 1:1-style ortholog rather than just the closest homolog)
# Caveat: a few gene models merge two neighbouring genes, so the top hit can
# describe only one of them.
BLASTP_OUT <- "/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.blastp.out"
GS_AT      <- "/private/groups/meyerlab/esteiger/inversions/synteny/genespace/syntenicHits/Deinandra_vs_Arabidopsis_thaliana.allBlast.txt.gz"
block_genes$mrna_id <- sub(".*mRNA_id=([^;]+).*", "\\1", genes$attributes[in_block])
bl <- read.table(BLASTP_OUT, sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
bl <- bl[bl$V1 %in% block_genes$mrna_id, ]
bl <- bl[order(bl$V1, -bl$V12), ]; bl <- bl[!duplicated(bl$V1), ]
block_genes$best_blast_hit <- bl$V2[match(block_genes$mrna_id, bl$V1)]
block_genes$best_blast_pident <- bl$V3[match(block_genes$mrna_id, bl$V1)]
block_genes$best_blast_description <- sub("^\\S+\\s+", "", bl$V14[match(block_genes$mrna_id, bl$V1)])
at <- read.table(gzfile(GS_AT), sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
at <- at[at$id1 %in% block_genes$mrna_id, ]
at <- at[order(at$id1, -at$bitScore), ]; at <- at[!duplicated(at$id1), ]
block_genes$arabidopsis_best_hit <- at$id2[match(block_genes$mrna_id, at$id1)]
block_genes$arabidopsis_same_orthogroup <- at$sameOG[match(block_genes$mrna_id, at$id1)]
#  - eggnog_name / eggnog_description: the eggNOG-mapper (Viridiplantae) run that
#    feeds the combined GO annotation; fills in genes with no blastp hit
EGGNOG <- "/private/groups/meyerlab/esteiger/genome/annotation/gvtp_eggnog_Di2_viridiplantae.emapper.annotations"
eg_lines <- readLines(EGGNOG)
eg_lines <- eg_lines[!grepl("^##", eg_lines)]
eg <- read.table(text = eg_lines, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
names(eg)[1] <- "query"
eg <- eg[eg$query %in% block_genes$mrna_id, ]
m <- match(block_genes$mrna_id, eg$query)
block_genes$eggnog_name <- ifelse(is.na(m) | eg$Preferred_name[m] == "-", NA, eg$Preferred_name[m])
block_genes$eggnog_description <- ifelse(is.na(m) | eg$Description[m] == "-", NA, eg$Description[m])
write.table(block_genes, file.path(output_dir, paste0("genes_", file_label, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
print(paste0("Genes in the block: ", nrow(block_genes), " (",
             paste(sapply(BLOCK$name, function(n) paste0(sum(block_genes$part == n), " in ", n)), collapse = "; "),
             "); with GO annotation: ", sum(block_genes$has_GO)))
print(paste0("Median gene length: block ", median(block_genes$length), " bp vs genome ", median(gene_length), " bp"))
target_genes <- block_genes$gene_id[block_genes$has_GO]


# --- Universe, length bias and categories ---
common_genes <- intersect(all_genes, names(gene_length))
gene.vector <- setNames(as.integer(common_genes %in% target_genes), common_genes)
bias.data <- gene_length[common_genes]
print(paste0("Universe: ", length(common_genes), " GO-annotated genes; targets in universe: ", sum(gene.vector)))

gene2cat <- gene2GO[common_genes]
cat_counts <- table(unlist(gene2cat))
keep_cats <- names(cat_counts)[cat_counts >= MIN_CAT_SIZE]
gene2cat <- lapply(gene2cat, function(x) x[x %in% keep_cats])
print(paste0("Category size filter (>=", MIN_CAT_SIZE, " genes): ", length(cat_counts), " -> ", length(keep_cats)))

jpeg(file.path(output_dir, paste0("pwf_gene_length_", file_label, ".jpg")), width = 800, height = 600)
pwf <- nullp(gene.vector, bias.data = bias.data, plot.fit = TRUE)
dev.off()


# --- Enrichment ---
GO.wall <- goseq(pwf, gene2cat = gene2cat, use_genes_without_cat = FALSE)
m <- match(GO.wall$category, term_tab$id)
GO.wall$term <- term_tab$name[m]; GO.wall$ontology <- term_tab$ontology[m]
GO.wall$FDR <- p.adjust(GO.wall$over_represented_pvalue, method = "BH")
GO.wall$Bonferroni <- p.adjust(GO.wall$over_represented_pvalue, method = "bonferroni")
gene_by_cat <- split(rep(names(gene2cat), lengths(gene2cat)), unlist(gene2cat))
GO.wall$block_genes <- sapply(GO.wall$category, function(go_id) {
    g <- gene_by_cat[[go_id]]; paste(g[g %in% target_genes], collapse = ", ")
})
GO.wall <- GO.wall[GO.wall$numDEInCat > 0, ]          # only terms with >= 1 block gene are informative
write.table(GO.wall, file.path(output_dir, paste0("GO_all_tested_", file_label, "_goseq.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
print("Top terms by nominal p (any ontology):")
print(head(GO.wall[order(GO.wall$over_represented_pvalue),
                   c("category", "ontology", "term", "numDEInCat", "numInCat", "over_represented_pvalue", "FDR")], 15),
      row.names = FALSE)

for (GROUP in c("BP", "MF", "CC")) {
    sig <- GO.wall[GO.wall$ontology == GROUP & !is.na(GO.wall$FDR) & GO.wall$FDR < 0.05, ]
    outfile <- file.path(output_dir, paste0("GO_enrichment_results_", file_label, "_", GROUP, "_goseq.tsv"))
    reduced_outfile <- sub("_goseq.tsv$", "_goseq_reduced.tsv", outfile)
    plots <- file.path(output_dir, paste0(c("GO_treemap_", "GO_scatterplot_"), file_label, "_", GROUP, "_goseq.jpg"))
    if (nrow(sig) > 0) {
        write.table(sig, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
        print(paste0(GROUP, ": ", nrow(sig), " FDR-significant terms -> ", outfile))
    } else {
        if (file.exists(outfile)) file.remove(outfile)
        print(paste0(GROUP, ": no FDR-significant terms"))
    }
    if (nrow(sig) >= 3) {
        simMatrix <- calculateSimMatrix(sig$category, orgdb = "org.At.tair.db", ont = GROUP, method = "Wang")
        reducedTerms <- reduceSimMatrix(simMatrix, scores = setNames(-log10(sig$FDR), sig$category),
                                        threshold = 0.7, orgdb = "org.At.tair.db")
        write.table(sig[sig$category %in% reducedTerms$go[reducedTerms$parent == reducedTerms$go], ],
                    reduced_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
        jpeg(plots[1], width = 1200, height = 800); treemapPlot(reducedTerms); dev.off()
        jpeg(plots[2], width = 900, height = 700); scatterPlot(simMatrix, reducedTerms); dev.off()
    } else {
        stale <- c(reduced_outfile, plots); file.remove(stale[file.exists(stale)])
    }
}
