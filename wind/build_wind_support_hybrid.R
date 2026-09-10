## Directed wind-support graph built from the HYBRID wind climatology
## (WindNinja direction + raw HRRR speed; see build_hybrid_windrose.py) --
## companion to wind_connectivity.R's section 4, factored out because the
## hybrid climatology is already per-population (not per HRRR-grid-cell), so
## the nearest-cell lookup wind_connectivity.R does is unnecessary here:
## cell_id in hybrid_windrose_*.txt IS the population code.
##
## Motivation: PTGC1 validation showed WindNinja's terrain downscaling
## improves wind DIRECTION at an exposed cape station (median error 16 deg ->
## 10 deg) but does not recover HRRR's speed underestimate there (28% -> 32%
## low, i.e. no better) -- so the wind_support() kernel below, which needs
## both a direction (which sector) and a speed (transport strength), should
## take direction from WindNinja and speed from HRRR rather than both from
## either source alone.

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(geosphere)

WIND_DIR <- "/private/groups/meyerlab/esteiger/wind"

name_map <- c(
  "Point_Conception" = "conception", "CDFW" = "cdfw",
  "VAFB_Arguello"    = "arguello",   "Strauss" = "strauss",
  "VAFB_Lion_Head"   = "lionhead",   "Dangermond" = "Dangermond"
)

field <- read.delim("/private/groups/meyerlab/esteiger/field_data.txt", header = TRUE, sep = "\t")
centroids <- aggregate(cbind(Lat, Long) ~ Clusters, data = field, FUN = mean)
centroids$pop <- name_map[centroids$Clusters]
rownames(centroids) <- centroids$pop
pops <- rownames(centroids)

coords <- as.matrix(centroids[pops, c("Long", "Lat")])
rownames(coords) <- pops

sectors <- read.delim(file.path(WIND_DIR, "hrrr", "hybrid_windrose_sectors.txt"), header = TRUE, sep = "\t")

wind_support <- function(pop, bearing_deg) {
  sec <- sectors[sectors$cell_id == pop, ]
  ang <- abs((sec$sector_center_deg - bearing_deg + 180) %% 360 - 180)
  w <- pmax(0, cos(ang * pi / 180))
  if (sum(w) == 0) return(0)
  w <- w / sum(w)
  sum(w * sec$frequency * sec$mean_speed_ms)
}

n <- length(pops)
wind_support_asym <- matrix(0, n, n, dimnames = list(pops, pops))
for (i in pops) {
  for (j in pops) {
    if (i == j) next
    b <- bearing(coords[i, ], coords[j, ]) %% 360
    wind_support_asym[i, j] <- wind_support(i, b)
  }
}

out <- file.path(WIND_DIR, "wind_support_asymmetric_hybrid.txt")
write.table(round(wind_support_asym, 4), out, sep = "\t", quote = FALSE)
cat("Saved hybrid directed wind-support matrix ->", out, "\n")

cat("\nAsymmetry check -- |wind_support(i,j) - wind_support(j,i)|, upper triangle:\n")
asym_diff <- abs(wind_support_asym - t(wind_support_asym))
print(round(asym_diff[upper.tri(asym_diff)], 3))
