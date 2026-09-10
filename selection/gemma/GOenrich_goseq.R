#!/usr/bin/env Rscript --vanilla

# Gene-length-corrected GO enrichment using goseq (Wallenius approximation).
# Unlike GOenrich.R (topGO/Fisher, no length correction), this fits a
# probability weighting function on gene length before testing, which
# matters here because longer genes have more chances to contain an
# outlier SNP by chance alone.
#
# Usage: Rscript GOenrich_goseq.R <target_file> [output_dir]

if (getRversion() < "4.4.0") `%||%` <- function(x, y) if (is.null(x)) y else x

library(goseq)
library(dplyr)
library(rrvgo)
library(GO.db)

# Minimum number of genes (in the universe) a GO category must have to be
# tested at all — mirrors topGO's nodeSize pruning. Categories with very few
# genes are statistically fragile: a couple of coincidentally-linked genes
# can produce an extreme p-value with almost no real power behind it.
MIN_CAT_SIZE <- 3

args <- commandArgs(trailingOnly = TRUE)
target_file <- if (length(args) >= 1) args[1] else "/private/groups/meyerlab/esteiger/selection/gemma/fog/gene_list_fog_40000_3_woInversion.txt"
output_dir  <- if (length(args) >= 2) args[2] else dirname(target_file)

file_label <- sub("^gene_list_|^genes_", "", sub("\\.txt$", "", basename(target_file)))

GFF <- "/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"


# --- Load GO annotation ---
# Source is the combined table (genome/annotation/merge_go_sources.sh): the
# eggNOG-mapper run (convert_eggnog_to_geneid_go.py, ~51% gene coverage,
# orthology-group-based) unioned with harvest_blast_go.py's NCBI Gene-level
# lookup for the ~15,937 genes eggNOG assigned zero GO terms to (~8,880
# recovered). Not the old annotated_sequences_lvl0_enrich_geneid_go.tsv —
# that file is transcript-level and needed a features.tracking-based
# collapse to gene level (see git history if that logic is ever needed
# again). This file is already gene-level (built directly from the GFF's
# gene_id), so no collapse step is needed — collapsing it through the old
# transcript2gene table would actually corrupt gene ids that contain a dot
# (e.g. "MSTRG.14" would get wrongly truncated to "MSTRG" by the old
# dot-strip fallback).
go_annot_raw <- read.table("/private/groups/meyerlab/esteiger/genome/annotated_sequences_combined_geneid_go.tsv",
                           header=TRUE, sep="\t", stringsAsFactors=FALSE)
names(go_annot_raw) <- c("gene", "go")
go_annot_raw$gene <- trimws(go_annot_raw$gene)
go_annot_raw <- unique(go_annot_raw[, c("gene", "go")])
print(paste0("Loaded ", nrow(go_annot_raw), " gene-GO pairs across ",
             length(unique(go_annot_raw$gene)), " genes"))

# Remap obsolete GO ids to their current successor where one exists — see
# GSEA_GO.R's comment and build_go_ancestor_map.py for how the remap table
# was built (only unambiguous `replaced_by` chains, not the ambiguous
# `consider:` alternatives).
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

# Drop annotations with no plausible plant equivalent (e.g. "embryonic
# skeletal system development") BEFORE propagation and testing. BLAST-hit-
# derived GO transfer (harvest_blast_go.py) isn't restricted to plant source
# organisms, so a handful of gene annotations carry over animal-only terms.
# See build_nonplant_go_exclusion_list.py for how this list was built. In
# practice this removes <1% of tested categories (checked empirically), so
# it's a data-quality fix, not expected to change which terms clear FDR.
excl_file <- "/private/groups/meyerlab/esteiger/genome/go_terms_excluded_nonplant.txt"
if (file.exists(excl_file)) {
    excl_ids <- read.table(excl_file, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")[[1]]
    n_before <- nrow(go_annot_raw)
    go_annot_raw <- go_annot_raw[!(go_annot_raw$go %in% excl_ids), ]
    print(paste0("Non-plant-plausible term filter: ", n_before, " -> ", nrow(go_annot_raw), " direct gene-GO pairs"))
} else {
    print(paste0("WARNING: non-plant exclusion list not found at ", excl_file, " — skipping that filter"))
}


# --- Propagate annotations up the GO hierarchy (true path rule) ---
# topGO does this automatically when it builds a topGOdata object: a gene
# annotated to a term is implicitly annotated to every ancestor of that term.
# goseq does not do this for a custom gene2cat, so we replicate it by hand.
# Ancestors come from a fresh go-basic.obo (build_go_ancestor_map.py), NOT
# GO.db's own GOBPANCESTOR/GOCCANCESTOR/GOMFANCESTOR tables — those were
# found to lag behind the current ontology (missing a recently-added is_a
# edge; see build_nonplant_go_exclusion_list.py for how that was caught).
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

go_annot <- go_annot_prop %>%
  group_by(gene) %>%
  summarise(GO_terms = paste(unique(go), collapse = ",")) %>%
  ungroup()

gene2GO <- strsplit(go_annot$GO_terms, ",")
names(gene2GO) <- go_annot$gene
all_genes <- names(gene2GO)


# --- Read the full GFF once: gene lengths + a transcript_id -> gene_id map ---
gff_all <- read.table(GFF, sep="\t", header=FALSE, quote="", comment.char="#",
                      stringsAsFactors=FALSE,
                      col.names=c("scaffold","source","feature","start","end",
                                  "score","strand","frame","attributes"))

# Gene lengths: genomic span (end - start + 1) of each "gene" feature.
# Genomic span (not CDS length) is the right proxy for length bias here,
# since outlier SNPs are genomic positions that can fall anywhere in or
# near a gene body, not just in coding sequence.
gff_genes <- gff_all[gff_all$feature == "gene", ]
gff_genes$gene_id <- sub(".*gene_id=([^;]+).*", "\\1", gff_genes$attributes)
gene_length <- setNames(gff_genes$end - gff_genes$start + 1, gff_genes$gene_id)
gff_gene_ids <- unique(gff_genes$gene_id)

# transcript_id -> gene_id, from every GFF row carrying both attributes
tx_ids <- sub(".*transcript_id=([^;]+).*", "\\1", gff_all$attributes)
gene_ids_col <- sub(".*gene_id=([^;]+).*", "\\1", gff_all$attributes)
has_both <- grepl("transcript_id=", gff_all$attributes) & grepl("gene_id=", gff_all$attributes)
tx2gene_gff <- setNames(gene_ids_col[has_both], tx_ids[has_both])
tx2gene_gff <- tx2gene_gff[!duplicated(names(tx2gene_gff))]


# --- Load target gene list ---
target_genes <- unique(trimws(scan(target_file, what="")))
target_genes <- sub("\\.$", "", target_genes)

# Safety net: gene_list_*.txt is normally already resolved by
# resolve_gene_ids.R at creation time, but if this is ever called on a list
# that bypassed that step, catch stray transcript IDs here too (e.g.
# transcript "g14399.t1" belongs to gene "g4525", not "g14399" — a bare
# suffix/prefix match against the extracted id misses that).
resolve_one <- function(g) {
    cand <- c(g, paste0(g, ".t1"), paste0(g, ".1"))
    hit <- cand[cand %in% names(tx2gene_gff)]
    if (length(hit) > 0) return(unname(tx2gene_gff[hit[1]]))
    pref <- grep(paste0("^", gsub("\\.", "\\\\.", g), "\\."), names(tx2gene_gff), value = TRUE)
    if (length(pref) > 0) return(unname(tx2gene_gff[pref[1]]))
    NA_character_
}
not_a_gene_id <- !(target_genes %in% gff_gene_ids)
resolvable <- vapply(target_genes[not_a_gene_id], resolve_one, character(1))
n_resolved <- sum(!is.na(resolvable))
if (n_resolved > 0) {
    print(paste0(n_resolved, " target entries were transcript IDs (not gene IDs) — resolved to their real gene"))
    target_genes[not_a_gene_id][!is.na(resolvable)] <- resolvable[!is.na(resolvable)]
    target_genes <- unique(target_genes)
}

targets <- length(target_genes)
conserved <- length(target_genes[target_genes %in% all_genes])
print(paste0("Number of targets was ", targets, ", but number found in gene universe was ", conserved))
target_genes <- target_genes[target_genes %in% all_genes]


# --- Align universe to genes with known length ---
common_genes <- intersect(all_genes, names(gene_length))
missing_len <- setdiff(all_genes, names(gene_length))
if (length(missing_len) > 0) {
    print(paste0(length(missing_len), " of ", length(all_genes),
                 " annotated genes had no length in the GFF and were dropped from the universe"))
}

gene.vector <- as.integer(common_genes %in% target_genes)
names(gene.vector) <- common_genes
bias.data <- gene_length[common_genes]
print(paste0("Genes marked as interesting in gene.vector: ", sum(gene.vector == 1)))

gene2cat <- gene2GO[common_genes]

# --- Drop GO categories with fewer than MIN_CAT_SIZE genes in the universe ---
cat_counts <- table(unlist(gene2cat))
keep_cats <- names(cat_counts)[cat_counts >= MIN_CAT_SIZE]
gene2cat <- lapply(gene2cat, function(x) x[x %in% keep_cats])
print(paste0("Category size filter (>=", MIN_CAT_SIZE, " genes): ",
             length(cat_counts), " categories -> ", length(keep_cats), " kept"))


# --- Fit probability weighting function on gene length ---
pwf <- nullp(gene.vector, bias.data = bias.data, plot.fit = FALSE)


# --- Run GO enrichment for BP, MF, and CC ---
GO.wall <- goseq(pwf, gene2cat = gene2cat, test.cats = c("GO:BP", "GO:MF", "GO:CC"))
GO.wall$FDR <- p.adjust(GO.wall$over_represented_pvalue, method = "BH")
GO.wall$Bonferroni <- p.adjust(GO.wall$over_represented_pvalue, method = "bonferroni")

gene_by_cat <- split(rep(names(gene2cat), lengths(gene2cat)), unlist(gene2cat))

for (GROUP in c("BP", "MF", "CC")) {
    ont_label <- switch(GROUP, BP = "BP", MF = "MF", CC = "CC")
    sub_res <- GO.wall[GO.wall$ontology == ont_label, ]

    sig_terms <- sub_res[!is.na(sub_res$FDR) & sub_res$FDR < 0.05, ]
    sig_terms$Genes <- sapply(sig_terms$category, function(go_id) {
        genes <- gene_by_cat[[go_id]]
        sig_genes <- genes[genes %in% target_genes]
        paste(sig_genes, collapse = ", ")
    })

    outfile <- file.path(output_dir, paste0("GO_enrichment_results_", file_label, "_", GROUP, "_goseq.tsv"))
    if (nrow(sig_terms) > 0) {
        write.table(sig_terms, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
        print(paste0("Wrote ", nrow(sig_terms), " significant terms to ", outfile))
    } else {
        # Remove any stale file from a previous run over the same label (e.g.
        # a prior method version, or a different candidate gene list) so a
        # null result here can't be mistaken for a leftover positive one
        if (file.exists(outfile)) file.remove(outfile)
        print(paste0("No FDR-significant terms for ", GROUP, " — not writing ", outfile))
    }

    # --- rrvgo: reduce redundancy among FDR-significant terms ---
    if (nrow(sig_terms) >= 3) {
        # method="Wang" is graph-topology-only and needs no annotation-frequency
        # (Information Content) data, unlike Resnik/Lin/Rel/Jiang — appropriate
        # here since org.At.tair.db is only a stand-in orgdb (no real orgdb
        # exists for this species) and we don't want to import Arabidopsis-
        # specific term-frequency bias into which terms look "informative"
        simMatrix <- calculateSimMatrix(sig_terms$category,
                                        orgdb = "org.At.tair.db",
                                        ont = GROUP,
                                        method = "Wang")
        reducedTerms <- reduceSimMatrix(simMatrix,
                                        scores = setNames(-log10(sig_terms$FDR), sig_terms$category),
                                        threshold = 0.7,
                                        orgdb = "org.At.tair.db")
        reduced_GO   <- reducedTerms[reducedTerms$parent == reducedTerms$go, ]
        final_results <- sig_terms[sig_terms$category %in% reduced_GO$go, ]

        reduced_outfile <- file.path(output_dir, paste0("GO_enrichment_results_", file_label, "_", GROUP, "_goseq_reduced.tsv"))
        write.table(final_results, file = reduced_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
        print(paste0("Wrote ", nrow(final_results), " reduced terms to ", reduced_outfile))

        jpeg(file.path(output_dir, paste0("GO_treemap_", file_label, "_", GROUP, "_goseq.jpg")), width=1200, height=800)
        treemapPlot(reducedTerms)
        dev.off()

        jpeg(file.path(output_dir, paste0("GO_scatterplot_", file_label, "_", GROUP, "_goseq.jpg")), width=900, height=700)
        scatterPlot(simMatrix, reducedTerms)
        dev.off()
    } else {
        # Clean up any stale reduced/plot files from a previous run that did
        # have >=3 significant terms, so a rerun that no longer does isn't
        # left with orphaned, no-longer-supported outputs
        stale <- c(
            file.path(output_dir, paste0("GO_enrichment_results_", file_label, "_", GROUP, "_goseq_reduced.tsv")),
            file.path(output_dir, paste0("GO_treemap_", file_label, "_", GROUP, "_goseq.jpg")),
            file.path(output_dir, paste0("GO_scatterplot_", file_label, "_", GROUP, "_goseq.jpg"))
        )
        file.remove(stale[file.exists(stale)])
        print(paste0("Skipping rrvgo for ", GROUP, ": only ", nrow(sig_terms), " FDR-significant terms (need >= 3)"))
    }
}
