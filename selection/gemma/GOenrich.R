#!/usr/bin/env Rscript --vanilla

# Define %||% for R < 4.4.0 (needed by packages compiled under R 4.4.1)
if (getRversion() < "4.4.0") `%||%` <- function(x, y) if (is.null(x)) y else x

library(topGO)
library(dplyr)
library(rrvgo)

# --- Configuration: accepts command-line arguments or falls back to defaults ---
# Usage: Rscript GOenrich.R <target_file> [output_dir]
args <- commandArgs(trailingOnly = TRUE)
target_file <- if (length(args) >= 1) args[1] else "/private/groups/meyerlab/esteiger/selection/gemma/4x/fog/gene_list_fog.txt"
output_dir  <- if (length(args) >= 2) args[2] else dirname(target_file)

# Derive a label from the filename for output files
file_label <- sub("^gene_list_|^genes_", "", sub("\\.txt$", "", basename(target_file)))


# --- Load gene-transcript mapping from features.tracking ---
# This maps transcript-level IDs to gene-level IDs so each gene is counted once
tracking <- read.table("/private/groups/meyerlab/esteiger/genome/features.tracking",
                       header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
transcript2gene <- setNames(trimws(tracking$Gene), trimws(tracking$Transcript))


# --- Load GO annotation (transcript-level IDs) ---
go_annot_raw <- read.table("/private/groups/meyerlab/esteiger/genome/annotated_sequences_lvl0_enrich_geneid_go.tsv",
                           header=TRUE, sep="\t", stringsAsFactors=FALSE)
names(go_annot_raw) <- c("transcript", "go")
go_annot_raw$transcript <- trimws(go_annot_raw$transcript)

# Map transcripts to gene-level IDs using tracking file;
# fall back to dot-stripping for any transcripts not in tracking
in_tracking <- go_annot_raw$transcript %in% names(transcript2gene)
go_annot_raw$gene <- go_annot_raw$transcript
go_annot_raw$gene[in_tracking] <- transcript2gene[go_annot_raw$transcript[in_tracking]]
go_annot_raw$gene[!in_tracking] <- trimws(sub("\\..*", "", go_annot_raw$transcript[!in_tracking]))

# Deduplicate: each gene counted once per GO term (even if it has multiple transcripts)
go_annot_raw <- unique(go_annot_raw[, c("gene", "go")])
print(paste0("Mapped ", sum(in_tracking), " of ", length(in_tracking),
             " annotation rows via tracking; ", sum(!in_tracking), " fell back to dot-strip"))


# Convert to topGO format (all GO terms per gene, comma-separated)
go_annot <- go_annot_raw %>%
  group_by(gene) %>%
  summarise(GO_terms = paste(unique(go), collapse = ",")) %>%
  ungroup()

gene2GO <- strsplit(go_annot$GO_terms, ",")
names(gene2GO) <- go_annot$gene
all_genes <- names(gene2GO)


# --- Load target gene list ---
target_genes <- unique(trimws(scan(target_file, what="")))
# Strip trailing dots (some gene lists have them)
target_genes <- sub("\\.$", "", target_genes)

# Report how many targets are in the gene universe
targets <- length(target_genes)
conserved <- length(target_genes[target_genes %in% all_genes])
print(paste0("Number of targets was ", targets, ", but number found in gene universe was ", conserved))
target_genes <- target_genes[target_genes %in% all_genes]


# Make factor for topGO (1 = interesting)
geneList <- factor(as.integer(all_genes %in% target_genes), levels = c(0, 1))
names(geneList) <- all_genes
print(paste0("Genes marked as interesting in geneList: ", sum(geneList == 1)))


# Run GO enrichment for BP, MF, and CC
for (GROUP in c("BP", "MF", "CC")) {
	GOdata <- new("topGOdata",
              ontology = GROUP,
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)

	resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

	allRes <- GenTable(GOdata, Fisher = resultFisher, orderBy = "Fisher", topNodes = 300)
	allRes$pvalue <- as.numeric(gsub("<", "", allRes$Fisher))
	# Apply FDR correction (Benjamini-Hochberg)
	allRes$FDR <- p.adjust(allRes$pvalue, method = "BH")
	# Apply Bonferroni correction
	allRes$Bonferroni <- p.adjust(allRes$pvalue, method = "bonferroni")

	# Filter for FDR-significant terms
	sig_terms <- allRes[allRes$FDR < 0.05, ]
	# Add gene names for each GO term
	sig_terms$Genes <- sapply(sig_terms$GO.ID, function(go_id) {
  	genes <- genesInTerm(GOdata, go_id)[[1]]
  	sig_genes <- genes[genes %in% target_genes]
  	paste(sig_genes, collapse = ", ")
	})

	outfile <- file.path(output_dir, paste0("GO_enrichment_results_", file_label, "_", GROUP, "_unfiltered.tsv"))
	write.table(sig_terms, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
	print(paste0("Wrote ", nrow(sig_terms), " significant terms to ", outfile))

	# --- rrvgo: reduce redundancy among FDR-significant terms ---
	fdr_terms <- allRes[allRes$FDR < 0.05, ]
	if (nrow(fdr_terms) >= 3) {
		simMatrix <- calculateSimMatrix(fdr_terms$GO.ID,
		                                orgdb = "org.At.tair.db",
		                                ont = GROUP,
		                                method = "Rel")
		reducedTerms <- reduceSimMatrix(simMatrix,
		                                scores = setNames(-log10(fdr_terms$FDR),
		                                                  fdr_terms$GO.ID),
		                                threshold = 0.7,
		                                orgdb = "org.At.tair.db")
		reduced_GO   <- reducedTerms[reducedTerms$parent == reducedTerms$go, ]
		final_results <- fdr_terms[fdr_terms$GO.ID %in% reduced_GO$go, ]

		reduced_outfile <- file.path(output_dir, paste0("GO_enrichment_results_", file_label, "_", GROUP, "_reduced.tsv"))
		write.table(final_results, file=reduced_outfile, sep="\t", quote=FALSE, row.names=FALSE)
		print(paste0("Wrote ", nrow(final_results), " reduced terms to ", reduced_outfile))

		jpeg(file.path(output_dir, paste0("GO_treemap_", file_label, "_", GROUP, ".jpg")), width=1200, height=800)
		treemapPlot(reducedTerms)
		dev.off()

		jpeg(file.path(output_dir, paste0("GO_scatterplot_", file_label, "_", GROUP, ".jpg")), width=900, height=700)
		scatterPlot(simMatrix, reducedTerms)
		dev.off()
	} else {
		print(paste0("Skipping rrvgo for ", GROUP, ": only ", nrow(fdr_terms), " FDR-significant terms (need >= 2)"))
	}
}
