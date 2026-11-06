#!/usr/bin/env Rscript
# Long-term linear trend in summer fog at each plant's location, from the
# Kling fog record (env/fog/Kling/timeseries.tif: one band per year,
# 1950-2014; the same file and extraction as the Kling_detrend code in
# env/extract_env.R). For each plant, a least-squares line is fitted to its
# annual Kling fog values ~ year; the slope measures how much fog has changed
# at that location over 1950-2014 (negative = fog loss).
#
# Output (env_var.txt is NOT modified): OUT, one row per SampleID with
#   kling_fog_trend_slope   change in Kling fog per year (Kling units)
#   kling_fog_trend_decade  the same per decade
#   kling_fog_trend_p       p-value of the slope
#   kling_fog_mean          mean of the annual values, 1950-2014
#
# Usage: Rscript fog_trend.R

.libPaths(c("/private/groups/meyerlab/esteiger/software/R_libs", .libPaths()))
suppressMessages(library(terra))

ENV_FILE <- "/private/groups/meyerlab/esteiger/env_var.txt"
KLING    <- "/private/groups/meyerlab/esteiger/env/fog/Kling/timeseries.tif"
OUT      <- "/private/groups/meyerlab/esteiger/env/fog_trend_by_sample.tsv"

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
