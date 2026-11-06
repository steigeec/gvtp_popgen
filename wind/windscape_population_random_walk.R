## Sensitivity check for the population tests in windscape_tests.R, using
## windscape random-walk deposition (pairwise_random_walk) instead of least-cost travel time.
##
## Inputs: windscape/hrrr_wind_rose.tif and hrrr_land_fraction.tif
##   (windscape_connectivity.R); pairwise population Fst; field_data.txt;
##   env_var.txt; population heterozygosity
## Outputs (in BASE_DIR/wind/windscape/random_walk/):
##   pop_deposition_<clim>_hl<h>.txt (row = source population)
##   windscape_random_walk_results.txt
## R libraries: R_libs_windscape, R_libs (windscape, terra, geosphere)
## Usage: Rscript windscape_population_random_walk.R

BASE_DIR <- "/private/groups/meyerlab/esteiger"
.libPaths(c(file.path(BASE_DIR, "software/R_libs_windscape"),
            file.path(BASE_DIR, "software/R_libs/")))
library(windscape)
library(terra)
library(geosphere)

WIND_DIR <- file.path(BASE_DIR, "wind")
WS_DIR <- file.path(WIND_DIR, "windscape")
RW_DIR <- file.path(WS_DIR, "random_walk")
dir.create(RW_DIR, showWarnings = FALSE)
HALF_LIVES <- c(0.5, 3, 24) # scenario values (h), not measured
WATER_WEIGHT <- 0.1         # conductance over water down-weighted by 90%
NPERM <- 9999

name_map <- c(
  "Point_Conception" = "conception", "CDFW" = "cdfw",
  "VAFB_Arguello"    = "arguello",   "Strauss" = "strauss",
  "VAFB_Lion_Head"   = "lionhead",   "Dangermond" = "Dangermond"
)

# ---------- Population data (as in windscape_tests.R) ----------

fst_df <- read.delim(
  file.path(BASE_DIR, "structure/Fst/populations/pairwise_auto/pairwise_Fst_populations.txt"),
  header = TRUE, sep = "\t")
pops <- sort(unique(c(fst_df$pop1, fst_df$pop2)))
fst_mat <- matrix(0, length(pops), length(pops), dimnames = list(pops, pops))
fst_mat[cbind(fst_df$pop1, fst_df$pop2)] <- fst_df$Fst_weighted
fst_mat[cbind(fst_df$pop2, fst_df$pop1)] <- fst_df$Fst_weighted
fst_lin <- fst_mat / (1 - fst_mat)

# The 13 excluded sample IDs have untrustworthy coordinates.
BAD_COORD_IDS <- c("4548", "4549", "4550", "4551", "4552", "4553", "4554",
                    "4555", "4556", "4557", "4559", "4561", "4562")
field <- read.delim(file.path(BASE_DIR, "field_data.txt"), header = TRUE, sep = "\t")
field$SampleID <- as.character(field$SampleID)
field <- field[!field$SampleID %in% BAD_COORD_IDS, ]
centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field, FUN = mean)
rownames(centroids) <- name_map[centroids$Clusters]
coords <- as.matrix(centroids[pops, c("Long", "Lat")])
geo_log <- log(distm(coords, fun = distHaversine) / 1000)
dimnames(geo_log) <- list(pops, pops)
diag(geo_log) <- 0

env_raw <- read.delim(file.path(BASE_DIR, "env_var.txt"), header = TRUE)
env_raw$SampleID <- as.character(env_raw$SampleID)
env_m <- merge(field[, c("SampleID", "Clusters")], env_raw[, c("SampleID", "elevation", "BIO18", "fog")], by = "SampleID")
for (v in c("elevation", "BIO18", "fog")) env_m[[v]][env_m[[v]] == -999999] <- NA
env_pop <- aggregate(cbind(elevation, BIO18, fog) ~ Clusters, data = env_m, FUN = function(x) mean(x, na.rm = TRUE))
rownames(env_pop) <- name_map[env_pop$Clusters]
env_dist <- as.matrix(dist(scale(env_pop[pops, c("elevation", "BIO18", "fog")])))

pop_div <- read.csv(file.path(BASE_DIR, "diversity/He/population/population_heterozygosity.csv"))
pop_div <- setNames(pop_div$heterozygosity, pop_div$population)[pops]


# ---------- Random walks and tests ----------

results <- list()
# Dangermond is excluded from the primary set (Fst inflated by low coverage).
pop_sets <- list(five_no_Dangermond = setdiff(pops, "Dangermond"), six = pops)
for (clim in "hrrr") {
  rose <- wind_rose(file.path(WS_DIR, sprintf("%s_wind_rose.tif", clim)))
  frac <- rast(file.path(WS_DIR, sprintf("%s_land_fraction.tif", clim)))
  rose <- weight_conductance(rose, ifel(frac >= 0.5, 1, WATER_WEIGHT))
  check_cell_distance(rose, coords)

  for (hl in HALF_LIVES) {
    dep <- pairwise_random_walk(rose, coords, half_life = hl)
    dimnames(dep) <- list(pops, pops)
    write.table(signif(dep, 6), file.path(RW_DIR, sprintf("pop_deposition_%s_hl%s.txt", clim, hl)),
                sep = "\t", quote = FALSE)
    off <- dep[row(dep) != col(dep)]
    cat(sprintf("%s, half-life %s h: off-diagonal deposition %.3g to %.3g (zeros: %d)\n",
                clim, hl, min(off), max(off), sum(off <= 0)))
    stopifnot(all(off > 0))

    # Wind distance = -ln(mean deposition), so a positive r supports wind.
    wind_dist <- -log(pairwise_means(dep))
    diag(wind_dist) <- 0
    for (set in names(pop_sets)) {
      s <- pop_sets[[set]]
      set.seed(42)
      # In div, [i, j] = ln(deposition j->i / deposition i->j).
      tests <- list(
        geo  = mantel_test(fst_lin[s, s], wind_dist[s, s], list(geo_log[s, s]), nperm = NPERM),
        env  = mantel_test(fst_lin[s, s], wind_dist[s, s], list(env_dist[s, s]), nperm = NPERM),
        both = mantel_test(fst_lin[s, s], wind_dist[s, s], list(geo_log[s, s], env_dist[s, s]), nperm = NPERM),
        div  = mantel_test(pairwise_ratios(pop_div[s]), t(pairwise_ratios(dep[s, s])), nperm = NPERM)
      )
      results[[length(results) + 1]] <- data.frame(
        climatology = clim, populations = set, half_life_h = hl,
        r_geo = tests$geo$stat, p_geo = tests$geo$p.value,
        r_env = tests$env$stat, p_env = tests$env$p.value,
        r_geo_env = tests$both$stat, p_geo_env = tests$both$p.value,
        r_diversity = tests$div$stat, p_diversity = tests$div$p.value)
    }
  }
}

res <- do.call(rbind, results)
write.table(res, file.path(RW_DIR, "windscape_random_walk_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, digits = 3)
