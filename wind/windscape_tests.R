## Isolation and diversity tests of wind connectivity with windscape (Kling &
## Ackerly 2021), at the population level (HRRR) and among Tranquillon Mtn
## individuals (100 m WindNinja).
##
## Inputs: pop_hours_hrrr.txt, strauss_hours_windninja.txt (windscape_connectivity.R);
##   population and individual pairwise Fst; population nucleotide diversity;
##   individual heterozygosity (het.txt); env_var.txt; field_data.txt
## Output: BASE_DIR/wind/windscape/windscape_test_results.txt
## R libraries: windscape 2.0.0 and terra 1.9-50 (CRAN) in BASE_DIR/software/R_libs_windscape,
##   installed with install.packages(c("windscape", "terra"), lib = <that dir>);
##   geosphere from BASE_DIR/software/R_libs
## Usage: Rscript windscape_tests.R
##
## Wind distance is -ln(mean of the two directional flows), so a positive r
## supports isolation by wind. Two-sided tests; 99,999 permutations for the
## Tranquillon isolation tests, 9,999 otherwise.
## Individual diversity is Asym + residual from the coverage model in
## het_plot.R, so that diversity ratios are defined.
## Pairs with identical coordinates have travel time floored at the smallest
## positive value (geographic distance floored at 1 m).

BASE_DIR <- "/private/groups/meyerlab/esteiger"
.libPaths(c(file.path(BASE_DIR, "software/R_libs_windscape"),
            file.path(BASE_DIR, "software/R_libs/")))
library(windscape)
library(geosphere)

WIND_DIR <- file.path(BASE_DIR, "wind")
WS_DIR <- file.path(WIND_DIR, "windscape")
NPERM <- 9999
NPERM_TRANQ <- 99999

name_map <- c(
  "Point_Conception" = "conception", "CDFW" = "cdfw",
  "VAFB_Arguello"    = "arguello",   "Strauss" = "strauss",
  "VAFB_Lion_Head"   = "lionhead",   "Dangermond" = "Dangermond"
)

read_hours <- function(file) {
  h <- as.matrix(read.delim(file.path(WS_DIR, file), row.names = 1, check.names = FALSE))
  off <- row(h) != col(h)
  h[off & h == 0] <- min(h[off & h > 0])
  h
}

env_raw <- read.delim(file.path(BASE_DIR, "env_var.txt"), header = TRUE)
env_raw$SampleID <- as.character(env_raw$SampleID)
field <- read.delim(file.path(BASE_DIR, "field_data.txt"), header = TRUE, sep = "\t")
field$SampleID <- as.character(field$SampleID)

fmt <- function(m) sprintf("r = %7.4f, p = %.4f", m$stat, m$p.value)

isolation_tests <- function(fst_lin, hours, geo_log, env_dist, nperm = NPERM) {
  wind_dist <- -log(pairwise_means(1 / hours))
  diag(wind_dist) <- 0
  set.seed(42)
  list(
    `Mantel Fst ~ geo dist`                      = mantel_test(fst_lin, geo_log, nperm = nperm),
    `Mantel Fst ~ wind dist`                     = mantel_test(fst_lin, wind_dist, nperm = nperm),
    `Partial Mantel Fst ~ wind | geo`            = mantel_test(fst_lin, wind_dist, list(geo_log), nperm = nperm),
    `Partial Mantel Fst ~ wind | environment`    = mantel_test(fst_lin, wind_dist, list(env_dist), nperm = nperm),
    `Partial Mantel Fst ~ wind | geo + environment` = mantel_test(fst_lin, wind_dist, list(geo_log, env_dist), nperm = nperm)
  )
}

diversity_test <- function(diversity, hours) {
  inflow_ratio <- t(pairwise_ratios(1 / hours))  # [i, j] = ln(flow j->i / flow i->j)
  div_ratio <- pairwise_ratios(diversity)        # [i, j] = ln(div_i / div_j)
  set.seed(42)
  list(`Diversity ratio ~ inflow/outflow ratio` = mantel_test(div_ratio, inflow_ratio, nperm = NPERM))
}

results <- list()

# ---------- 1. Populations ----------

fst_df <- read.delim(
  file.path(BASE_DIR, "structure/Fst/populations/pairwise_auto/pairwise_Fst_populations.txt"),
  header = TRUE, sep = "\t"
)
pops <- sort(unique(c(fst_df$pop1, fst_df$pop2)))
n <- length(pops)
stopifnot(nrow(fst_df) == n * (n - 1) / 2)
fst_mat <- matrix(0, n, n, dimnames = list(pops, pops))
fst_mat[cbind(fst_df$pop1, fst_df$pop2)] <- fst_df$Fst_weighted
fst_mat[cbind(fst_df$pop2, fst_df$pop1)] <- fst_df$Fst_weighted
pop_fst_lin <- fst_mat / (1 - fst_mat)  # Rousset 1997 linearization

# 13 individuals with untrustworthy coordinates are excluded from centroids.
BAD_COORD_IDS <- c("4548", "4549", "4550", "4551", "4552", "4553", "4554",
                    "4555", "4556", "4557", "4559", "4561", "4562")
field_ok <- field[!field$SampleID %in% BAD_COORD_IDS, ]
centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field_ok, FUN = mean)
rownames(centroids) <- name_map[centroids$Clusters]
pop_geo_log <- log(distm(as.matrix(centroids[pops, c("Long", "Lat")]), fun = distHaversine) / 1000)
diag(pop_geo_log) <- 0
dimnames(pop_geo_log) <- list(pops, pops)

# -999999 is env_var.txt's missing-data sentinel.
env_m <- merge(field_ok[, c("SampleID", "Clusters")], env_raw[, c("SampleID", "elevation", "BIO18", "fog")], by = "SampleID")
for (v in c("elevation", "BIO18", "fog")) env_m[[v]][env_m[[v]] == -999999] <- NA
env_pop <- aggregate(cbind(elevation, BIO18, fog) ~ Clusters, data = env_m, FUN = function(x) mean(x, na.rm = TRUE))
rownames(env_pop) <- name_map[env_pop$Clusters]
pop_env_dist <- as.matrix(dist(scale(env_pop[pops, c("elevation", "BIO18", "fog")])))

pop_div <- read.csv(file.path(BASE_DIR, "diversity/He/population/population_heterozygosity.csv"))
pop_div <- setNames(pop_div$heterozygosity, pop_div$population)[pops]
stopifnot(!anyNA(pop_div))

pop_hours <- read_hours("pop_hours_hrrr.txt")
pop_sets <- list(`5, Dangermond excluded` = setdiff(pops, "Dangermond"), `6` = pops)
for (set in names(pop_sets)) {
  s <- pop_sets[[set]]
  results[[sprintf("Populations (n = %s), HRRR wind", set)]] <-
    c(isolation_tests(pop_fst_lin[s, s], pop_hours[s, s], pop_geo_log[s, s], pop_env_dist[s, s]),
      diversity_test(pop_div[s], pop_hours[s, s]))
}

# ---------- 2. Tranquillon Mtn individuals, WindNinja wind ----------

strauss_hours <- read_hours("strauss_hours_windninja.txt")

# Isolation: individuals with pairwise Fst.
bam_ids <- sub("_rmdup\\.bam$", "", readLines(file.path(BASE_DIR, "mapping/merged/bamlist_strauss.txt")))
fst_all <- read.delim(file.path(BASE_DIR, "structure/Fst/individuals/pairwise_Fst.txt"), header = TRUE)
fst_all$sample1 <- as.character(fst_all$sample1); fst_all$sample2 <- as.character(fst_all$sample2)
fst_s <- fst_all[fst_all$sample1 %in% bam_ids & fst_all$sample2 %in% bam_ids, ]
# 6432 duplicates 5056 (same coordinates, extreme Fst); dropped.
fst_s <- fst_s[fst_s$sample1 != "6432" & fst_s$sample2 != "6432", ]
ids <- sort(unique(c(fst_s$sample1, fst_s$sample2)))
stopifnot(nrow(fst_s) == length(ids) * (length(ids) - 1) / 2, all(ids %in% rownames(strauss_hours)))
fst_mat <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
fst_mat[cbind(fst_s$sample1, fst_s$sample2)] <- fst_s$Fst_weighted
fst_mat[cbind(fst_s$sample2, fst_s$sample1)] <- fst_s$Fst_weighted
ind_fst_lin <- fst_mat / (1 - fst_mat)

geo_m <- distm(as.matrix(field[match(ids, field$SampleID), c("Long", "Lat")]), fun = distHaversine)
# One pair shares identical GPS coordinates; floor at 1 m before logging.
geo_m[geo_m == 0 & row(geo_m) != col(geo_m)] <- 1
ind_geo_log <- log(geo_m / 1000)
diag(ind_geo_log) <- 0
env_ind <- env_raw[match(ids, env_raw$SampleID), c("elevation", "BIO18", "fog")]
stopifnot(!anyNA(env_ind))  # complete for all Strauss individuals
ind_env_dist <- as.matrix(dist(scale(env_ind)))

# Diversity: coverage-corrected individual heterozygosity (see header).
het <- read.table(file.path(BASE_DIR, "diversity/He/individual/het.txt"), col.names = c("Sample", "Het"))
het$Sample <- sub("_rmdup\\.ml$", "", het$Sample)
het <- merge(het, field[, c("SampleID", "Cov")], by.x = "Sample", by.y = "SampleID")
het <- het[het$Sample != "4531", ]  # excluded in het_plot.R (likely contaminated)
cov_mod <- nls(Het ~ SSasymp(Cov, Asym, R0, lrc), data = het)
het$corrected <- coef(cov_mod)[["Asym"]] + residuals(cov_mod)
individuals <- read.delim(file.path(WIND_DIR, "strauss_individual", "strauss_individuals.txt"), header = TRUE)
div_ids <- as.character(individuals$id)
ind_div <- setNames(het$corrected, het$Sample)[div_ids]
stopifnot(!anyNA(ind_div), all(ind_div > 0))

results[["Tranquillon Mtn individuals, WINDNINJA wind"]] <- c(
  isolation_tests(ind_fst_lin, strauss_hours[ids, ids], ind_geo_log, ind_env_dist, NPERM_TRANQ),
  diversity_test(ind_div, strauss_hours[div_ids, div_ids])
)

# ---------- 3. Summary ----------

summary_txt <- file.path(WS_DIR, "windscape_test_results.txt")
sink(summary_txt)
cat("Wind-genetics hypothesis tests -- windscape (Kling & Ackerly 2021), two-sided, ", NPERM, " permutations (Tranquillon isolation: ", NPERM_TRANQ, ")\n", sep = "")
cat("ISOLATION: linearized Fst ~ -ln(mean wind flow); positive r = weaker wind connectivity, more differentiation\n")
cat("DIVERSITY: ln(div_i/div_j) ~ ln(flow j->i / flow i->j); positive r = downwind (net-receiving) sites more diverse\n")
cat(sprintf("Tranquillon isolation n = %d individuals; diversity n = %d individuals\n", length(ids), length(div_ids)))
cat("=====================================================================================================\n")
for (set in names(results)) {
  cat("\n", set, "\n", sep = "")
  for (test in names(results[[set]])) cat(sprintf("  %-48s %s\n", test, fmt(results[[set]][[test]])))
}
sink()
cat(readLines(summary_txt), sep = "\n")
cat("\nSaved summary ->", summary_txt, "\n")
