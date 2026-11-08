#!/usr/bin/env Rscript
# Linear trend in summer fog at each plant's location, from the annual Kling
# fog record (env/fog/Kling/timeseries.tif, 1950-2014). Fits annual fog ~ year
# per plant; negative slope = fog loss.
#
# Output: OUT, one row per SampleID (env_var.txt is not modified) with
#   kling_fog_trend_slope   change in Kling fog per year
#   kling_fog_trend_decade  the same per decade
#   kling_fog_trend_p       p-value of the slope
#   kling_fog_mean          mean of the annual values, 1950-2014
#
# Usage: Rscript fog_trend.R

BASE_DIR <- "/private/groups/meyerlab/esteiger"

.libPaths(c(file.path(BASE_DIR, "software/R_libs"), .libPaths()))
suppressMessages(library(terra))

ENV_FILE <- file.path(BASE_DIR, "env_var.txt")
KLING    <- file.path(BASE_DIR, "env/fog/Kling/timeseries.tif")
OUT      <- file.path(BASE_DIR, "env/fog_trend_by_sample.tsv")

env <- read.delim(ENV_FILE, colClasses = c(SampleID = "character"))
env <- env[!duplicated(env$SampleID) & !is.na(env$Lat) & !is.na(env$Long), ]
kr <- rast(KLING)
years <- as.integer(names(kr))
pts <- project(vect(env, geom = c("Long", "Lat"), crs = "EPSG:4326"), crs(kr))
annual <- as.matrix(extract(kr, pts)[, -1])              # plants x years
cat(sprintf("Kling record: %d years (%d-%d); units example: %s\n", length(years), min(years), max(years),
            paste(round(annual[1, 1:3], 3), collapse = ", ")))

fit <- t(apply(annual, 1, function(a) {
    ok <- !is.na(a)
    if (sum(ok) < 10) return(c(NA, NA))
    s <- summary(lm(a[ok] ~ years[ok]))$coefficients
    c(s[2, 1], s[2, 4])
}))
out <- data.frame(SampleID = env$SampleID, Location = env$Location,
                  kling_fog_trend_slope = fit[, 1], kling_fog_trend_decade = fit[, 1] * 10,
                  kling_fog_trend_p = fit[, 2], kling_fog_mean = rowMeans(annual, na.rm = TRUE))
write.table(out, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote %d plants to %s\n", nrow(out), OUT))
print(aggregate(cbind(kling_fog_mean, kling_fog_trend_decade, kling_fog_trend_p) ~ Location, data = out, FUN = mean),
      digits = 3)
