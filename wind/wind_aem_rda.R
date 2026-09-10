## AEM + RDA: directional test of wind connectivity on population allele
## frequencies (Blanchet et al. 2008, "Modelling directional spatial
## processes in ecological data").
##
## WHY THIS EXISTS (see structure/mantel_IBD.R and wind/wind_connectivity.R
## for the non-directional precursors): Fst is symmetric by construction, so
## no test built on pairwise Fst can ever show that wind-mediated dispersal
## is directional (e.g. that gene flow runs disproportionately downwind).
## AEM sidesteps that by changing the RESPONSE from pairwise Fst to the
## population allele-frequency matrix itself (populations x loci), and
## representing wind as a genuinely DIRECTED graph among populations (edge
## i->j gets its own weight, generally != edge j->i). AEM eigenvectors
## summarize that directed graph as a set of orthogonal directional spatial
## predictors, which are then tested against allele frequencies by RDA with
## geographic distance partialled out -- so a significant result means wind
## direction/asymmetry explains genetic variation beyond what geography
## alone would, not just "these sites are windy".
##
## This avoids the divMigrate/Gst-asymmetry alternative's known confound
## (asymmetric drift/Ne can mimic asymmetric migration in that method) --
## AEM only uses the (already symmetric-Fst-based, standard) allele
## frequencies as the response; directionality enters purely through the
## predictor side, where we supply it directly from real wind data instead
## of inferring it from the genetic data itself.
##
## DATA SOURCE: genotype-likelihood-based per-population allele frequencies
## (ANGSD -doMaf 1, reading GLs directly from the ANGSD-native BCF via
## -vcf-PL; see build_pop_allele_freq.sh Part 1), NOT the PLINK hard-called
## genotypes used previously. This fixes the earlier Dangermond problem:
## PLINK's 4x-coverage filter left only 1 of Dangermond's 4 individuals,
## making its "allele frequency" just one individual's genotype; the GL-based
## approach uses all 4 (matching the SAF/Fst pipeline's sample size) and
## additionally uses every individual regardless of coverage for every OTHER
## population too (down to 0.6x), not just the >4x subset PLINK required.
## Dangermond (n=4) is still the smallest, noisiest population, so both the
## 6-population and 5-population (Dangermond excluded) results are still
## reported -- but with real n=4 data rather than a single-genotype
## artifact, the 6-population version is now treated as primary, and the
## 5-population version is a robustness check rather than a required
## exclusion.
##
## Every population's ANGSD run used -doMajorMinor 4 (reference allele =
## major, matching -anc/-ref = same reference genome used in pop_SAF.sh),
## so "major" is trivially identical across all 6 populations for every site
## (it's just the reference base) -- no per-population major/minor
## reconciliation is needed; each file's "knownEM" column is already a
## directly comparable frequency of the non-reference allele.
##
## Not every site appears in every population's .mafs.gz: ANGSD's default
## SNP significance filter (-SNP_pval, lenient default) drops sites where a
## given population shows no evidence of polymorphism -- e.g. Dangermond
## (n=4, least power) retains 36,450,246 of 36,765,900 total sites, Strauss
## (n=738) retains nearly all of them. A site missing from a population's
## file means that population's own data couldn't reject "fixed for the
## reference allele" there, so missing values are filled with 0 (not NA)
## when building the combined matrix below -- the defensible reading of
## "not detected as polymorphic," and essential for not losing
## population-private variants (a site private to one small population
## would otherwise be dropped by an inner join, which is exactly the signal
## this analysis is trying to detect).

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(data.table)
library(vegan)
library(adespatial)
library(geosphere)

WIND_DIR <- "/private/groups/meyerlab/esteiger/wind"
# Raw/derived allele-frequency data lives with the other raw wind downloads
# in wind/, not the git-tracked project directory.
ALLELE_FREQ_DIR <- "/private/groups/meyerlab/esteiger/wind/allele_freq"

# ---------- 1. Population allele frequencies (GL-based, ANGSD -doMaf;
#               see build_pop_allele_freq.sh Part 1) ----------

maf_files <- c(
  CDFW             = "CDFW_maf.mafs.gz",
  Dangermond       = "Dangermond_maf.mafs.gz",
  Point_Conception = "Point_Conception_maf.mafs.gz",
  Strauss          = "Strauss_maf.mafs.gz",
  VAFB_Arguello    = "VAFB_Arguello_maf.mafs.gz",
  VAFB_Lion_Head   = "VAFB_Lion_Head_maf.mafs.gz"
)

t0 <- Sys.time()
pop_tables <- list()
for (pop in names(maf_files)) {
  dt <- fread(file.path(ALLELE_FREQ_DIR, maf_files[pop]),
              select = c("chromo", "position", "knownEM"))
  setnames(dt, "knownEM", pop)
  setkey(dt, chromo, position)
  pop_tables[[pop]] <- dt
  cat(sprintf("  %-17s %d sites [%.0fs]\n", pop, nrow(dt),
              as.numeric(Sys.time() - t0, units = "secs")))
}

# Full outer join across all 6 populations on (chromo, position); see header
# for why missing values are filled with 0 rather than left NA or inner-joined away.
wide <- Reduce(function(a, b) merge(a, b, by = c("chromo", "position"), all = TRUE), pop_tables)
for (pop in names(maf_files)) {
  set(wide, which(is.na(wide[[pop]])), pop, 0)
}
cat(sprintf("Merged wide table: %d sites x %d populations [%.0fs total]\n",
            nrow(wide), length(maf_files), as.numeric(Sys.time() - t0, units = "secs")))

mat <- as.matrix(wide[, names(maf_files), with = FALSE])
rownames(mat) <- paste(wide$chromo, wide$position, sep = ":")
Y_all <- t(mat)  # populations (rows) x loci (columns)
storage.mode(Y_all) <- "double"

# PLINK FID (= raw field_data.txt Clusters value) uses e.g. "VAFB_Arguello",
# "Point_Conception"; wind_support_asymmetric.txt (written by
# wind_connectivity.R) uses the short Fst-style names ("arguello",
# "conception", ...). Rename rows here to the short names so both files line
# up, matching the convention already used throughout structure/mantel_IBD.R
# and wind_connectivity.R.
name_map <- c(
  "Point_Conception" = "conception",
  "CDFW"             = "cdfw",
  "VAFB_Arguello"    = "arguello",
  "Strauss"          = "strauss",
  "VAFB_Lion_Head"   = "lionhead",
  "Dangermond"       = "Dangermond"
)
rownames(Y_all) <- name_map[rownames(Y_all)]

cat(sprintf("Allele frequency matrix: %d populations x %d loci\n", nrow(Y_all), ncol(Y_all)))

# MAF filter: keep a locus only if at least one population's own frequency
# exceeds MIN_MAX_MAF, rather than filtering on the POOLED frequency across
# all individuals (the existing postParalog_MAF_auto.sites approach, which
# is really "common in Strauss" since Strauss is 88% of the pooled sample --
# see chat discussion, deliberately not used here). This keeps
# population-private signal (e.g. a Dangermond-private allele at 20% passes
# easily, max=0.20) while dropping sites where every population sits in the
# noise floor -- important because sample sizes are wildly unequal
# (Dangermond n=4 = 8 allele copies -> any observed copy is already 12.5%,
# so this filter is a no-op there; Strauss n=738 = 1476 allele copies ->
# resolution of 0.07%, fine enough to "detect" huge numbers of singleton/
# doubleton variants that are sampling noise, not real differentiation, and
# that drop_monomorphic() alone would happily keep since var>0 is a very
# weak bar).
MIN_MAX_MAF <- 0.05

filter_sites <- function(Y, label) {
  max_freq <- apply(Y, 2, max)
  keep_maf <- max_freq > MIN_MAX_MAF
  cat(sprintf("  %s: %d / %d loci pass max-population-MAF > %.2f\n",
              label, sum(keep_maf), ncol(Y), MIN_MAX_MAF))
  Y <- Y[, keep_maf, drop = FALSE]
  keep_poly <- apply(Y, 2, function(x) var(x) > 0)
  cat(sprintf("  %s: %d / %d of those also polymorphic across these %d populations (kept)\n",
              label, sum(keep_poly), ncol(Y), nrow(Y)))
  Y[, keep_poly, drop = FALSE]
}

cat("Filtering loci (all 6 populations):\n")
Y_all <- filter_sites(Y_all, "6-pop")

Y_no_dang <- Y_all[rownames(Y_all) != "Dangermond", ]
cat("Filtering loci (5 populations, Dangermond excluded):\n")
Y_no_dang <- filter_sites(Y_no_dang, "5-pop")

# ---------- 2. Site centroids + geographic distance ----------

field <- read.delim("/private/groups/meyerlab/esteiger/field_data.txt", header = TRUE, sep = "\t")
centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field, FUN = mean)
centroids$pop <- name_map[centroids$Clusters]
rownames(centroids) <- centroids$pop

coords_all <- as.matrix(centroids[rownames(Y_all), c("Long", "Lat")])
geo_all <- as.matrix(distm(coords_all, fun = distHaversine)) / 1000
dimnames(geo_all) <- list(rownames(Y_all), rownames(Y_all))

# ---------- 3. Directed wind graph -> AEM eigenvectors ----------
# wind_support_asymmetric.txt was written by wind_connectivity.R: rows =
# origin population, columns = destination, entries = frequency/speed
# weighted downwind transport score FROM row TO column (not symmetric).

wind_asym_full <- as.matrix(read.delim(file.path(WIND_DIR, "wind_support_asymmetric.txt"), row.names = 1))

build_aem <- function(pops) {
  w <- wind_asym_full[pops, pops]
  n <- length(pops)
  # One directed edge per ordered pair (i != j): n*(n-1) edges total.
  # Site-by-edge incidence: +1 at the edge's origin, -1 at its destination.
  edges <- expand.grid(from = pops, to = pops, stringsAsFactors = FALSE)
  edges <- edges[edges$from != edges$to, ]
  edges$weight <- mapply(function(i, j) w[i, j], edges$from, edges$to)
  se_mat <- matrix(0, nrow = n, ncol = nrow(edges), dimnames = list(pops, NULL))
  for (e in seq_len(nrow(edges))) {
    se_mat[edges$from[e], e] <- 1
    se_mat[edges$to[e], e] <- -1
  }
  # Edge weight normalized to (0,1]: stronger downwind transport -> stronger
  # AEM connection weight, following the standard AEM weighting convention
  # (Blanchet et al. 2008 eq. for w = 1 - (d/dmax) adapted here to wind
  # support directly, since higher wind_support already means "more
  # connected" rather than "more distant").
  edge_w <- edges$weight / max(edges$weight)
  ae <- aem(binary.mat = se_mat, weight = edge_w)
  list(vectors = ae$vectors, values = ae$values, edges = edges, se_mat = se_mat)
}

cat("\n=== Building AEM eigenvectors (6 populations) ===\n")
aem_all <- build_aem(rownames(Y_all))
cat("AEM eigenvalues (6 pop):\n"); print(round(aem_all$values, 4))

cat("\n=== Building AEM eigenvectors (5 populations, no Dangermond) ===\n")
aem_nd <- build_aem(rownames(Y_no_dang))
cat("AEM eigenvalues (5 pop):\n"); print(round(aem_nd$values, 4))

# ---------- 4. RDA: allele frequencies ~ AEM1 | Condition(geo PCoA1) ----------
# n is tiny (5-6 populations) -- only the single strongest AEM axis is
# tested, conditioning on the first geographic PCoA axis as the "plain
# distance" null. This matches the total degrees of freedom actually
# available (n-1) rather than overfitting a model with more terms than the
# data can support.

run_rda <- function(Y, geo, aem_vectors, label) {
  geo_pcoa <- cmdscale(as.dist(geo), k = 1)
  aem1 <- aem_vectors[, 1, drop = FALSE]
  colnames(aem1) <- "AEM1"
  df <- data.frame(geoPCoA1 = geo_pcoa[, 1], AEM1 = aem1[, 1])
  rownames(df) <- rownames(Y)

  mod <- rda(Y ~ AEM1 + Condition(geoPCoA1), data = df)
  set.seed(42)
  t0 <- Sys.time()
  # n is tiny (5-6 pops), so vegan/permute will detect that the requested
  # count exceeds the exact permutation null and fall back to complete
  # enumeration (as happened automatically with the Mantel tests earlier,
  # "Set of permutations < minperm. Generating entire set."). Requesting a
  # moderate number (not 9999) avoids one failure mode: if the exact count
  # happens to exceed a very low request, vegan would silently do a random
  # subset instead of exact enumeration.
  test <- anova.cca(mod, permutations = 4999, by = NULL)
  cat(sprintf("[RDA + permutation test: %.0fs]\n", as.numeric(Sys.time() - t0, units = "secs")))
  cat(sprintf("\n--- RDA: %s ---\n", label))
  print(mod)
  cat("\nPermutation test (AEM1 | geoPCoA1):\n")
  print(test)
  list(model = mod, test = test)
}

cat("\n############ 6-population model (PRIMARY -- Dangermond now n=4, GL-based) ############\n")
res_all <- run_rda(Y_all, geo_all, aem_all$vectors, "6 populations")

cat("\n############ 5-population model (Dangermond excluded -- robustness check) ############\n")
geo_nd <- geo_all[rownames(Y_no_dang), rownames(Y_no_dang)]
res_nd <- run_rda(Y_no_dang, geo_nd, aem_nd$vectors, "5 populations, no Dangermond")

# ---------- 5. Save summary ----------

summary_txt <- file.path(WIND_DIR, "wind_aem_rda_results.txt")
sink(summary_txt)
cat("AEM + RDA directional wind connectivity test -- summary\n")
cat("=========================================================\n\n")
cat("PRIMARY result: 6 populations (GL-based allele frequencies; Dangermond n=4, see script header)\n")
print(res_all$test)
cat("\n\nRobustness check: 5 populations, Dangermond excluded (smallest/noisiest population)\n")
print(res_nd$test)
sink()
cat("\nSaved summary ->", summary_txt, "\n")
