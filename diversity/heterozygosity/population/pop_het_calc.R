# Population expected heterozygosity from folded SFS (.ml) files, and comparison
# with nucleotide diversity across plant taxa (Leffler et al. dataset).

BASE_DIR <- "/private/groups/meyerlab/esteiger"

setwd(file.path(BASE_DIR, "diversity/He/population"))

populations <- c("cdfw", "lionhead", "strauss", "arguello", "Dangermond", "conception")

results <- data.frame(population = character(), heterozygosity = numeric(),
                      stringsAsFactors = FALSE)

for (pop in populations) {
  filename <- paste0(pop, ".ml")

  if (!file.exists(filename)) {
    message("Skipping ", pop, ": ", filename, " not found")
    next
  }

  a <- scan(filename, quiet = TRUE)

  # The folded SFS has n_ind + 1 entries (R indices 1 .. n_ind+1):
  #   a[1]     = invariable sites (0 minor allele copies)
  #   a[k+1]   = sites with k minor allele copies, k = 1 .. n_ind
  # Total haplotypes = 2 * n_ind = 2 * (length(a) - 1)
  n_hap <- (length(a) - 1) * 2

  # Minor allele frequency for each polymorphic bin (k = 1 .. n_ind)
  freqs <- seq_len(length(a) - 1) / n_hap

  # Proportion of all sites (variable + invariable) in each polymorphic bin
  site_freqs <- a[-1] / sum(a)

  # Expected heterozygosity: sum of 2p(1-p) * f across frequency classes
  theta <- sum(2 * freqs * (1 - freqs) * site_freqs)

  results <- rbind(results, data.frame(population = pop, heterozygosity = theta,
                                       stringsAsFactors = FALSE))
  message("For ", pop, ", heterozygosity = ", round(theta, 6))
}

print(results)

## --- Compare nucleotide diversity across taxa (Leffler et al. plant dataset) ---

library(dplyr)
library(tidyr)
library(ggplot2)

leffler <- read.delim("plant_Leffler_nucleotide_div.txt", stringsAsFactors = FALSE)
colnames(leffler) <- c("latin_name", "common_name", "phylum", "class", "total_loci",
                        "diversity_pct", "theta_or_pi", "chrom_type", "site_type",
                        "sampling_strategy", "reference_pmid", "mating_system",
                        "habitat", "range_category", "included_main", "included_sex_auto")

# One dot per population: split the semicolon-separated per-population values
leffler_dots <- leffler %>%
  transmute(taxon = latin_name, phylum, mating_system, range_category,
            diversity = diversity_pct) %>%
  separate_rows(diversity, sep = ";") %>%
  mutate(diversity = as.numeric(diversity)) %>%
  filter(!is.na(diversity))

# Our study taxon: one dot per sampled population (SFS-derived heterozygosity, as %)
own_taxon <- "Deinandra increscens villosa"
own_dots <- data.frame(
  taxon = own_taxon,
  phylum = "Magnoliophyta",
  mating_system = "Obligatory outcrossing",
  range_category = "Narrow endemic",
  diversity = results$heterozygosity * 100,
  stringsAsFactors = FALSE
)

dot_data <- bind_rows(leffler_dots, own_dots)

# Sort taxa (y-axis) by their mean nucleotide diversity
taxon_order <- dot_data %>%
  group_by(taxon) %>%
  summarize(mean_diversity = mean(diversity), .groups = "drop") %>%
  arrange(mean_diversity)

dot_data <- dot_data %>%
  mutate(taxon = factor(taxon, levels = taxon_order$taxon))

# Panel per phylum (top to bottom): Angiospermae, Coniferae, Chlorophyta
dot_data <- dot_data %>%
  mutate(phylum_group = case_when(
           phylum == "Pinophyta"    ~ "Coniferae",
           phylum == "Chlorophyta"  ~ "Chlorophyta",
           TRUE                     ~ "Angiospermae"
         ),
         phylum_group = factor(phylum_group,
                                levels = c("Angiospermae", "Coniferae", "Chlorophyta")))

# Fill color by range category
range_colors <- c(
  "Narrow endemic" = "#FFFF00",
  "Broad endemic"  = "#4CAF50",
  "Cosmopolitan"   = "#145214",
  "Unknown"        = "grey60"
)

# Point shape by mating system: circle = outcrossing, triangle = partial selfing
mating_shapes <- c(
  "Obligatory outcrossing" = 21,
  "Partial selfing"        = 24
)

# Y-axis labels: italic Latin name only (bold italic for our taxon)
taxon_label_exprs <- lapply(as.character(taxon_order$taxon), function(nm) {
  if (nm == own_taxon) bquote(bolditalic(.(nm))) else bquote(italic(.(nm)))
})
names(taxon_label_exprs) <- taxon_order$taxon

p <- ggplot(dot_data, aes(x = diversity, y = taxon, fill = range_category, shape = mating_system)) +
  geom_point(color = "black", stroke = 0.4, size = 3.5) +
  scale_fill_manual(values = range_colors, name = "Range category") +
  scale_shape_manual(values = mating_shapes, name = "Mating system") +
  scale_y_discrete(labels = function(x) taxon_label_exprs[x]) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  facet_grid(phylum_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Nucleotide diversity (%)", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey85", color = NA),
    strip.text.y.left = element_text(size = 7, angle = 90)
  )

ggsave("nucleotide_diversity_dotplot.jpg", p, width = 9, height = 9, dpi = 300)
