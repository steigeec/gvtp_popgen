#!/usr/bin/env python3
"""
Download Global Wind Atlas mean wind speed (10 m height, 250 m resolution)
for the Point Conception study area, and extract site-level values at the
six population centroids.

Why GWA and not TerraClimate wind (already in env_var.txt / select_env_var.R):
TerraClimate wind speed is ~4.6 km resolution -- coarser than the ~5-20 km
spacing between several of our sites, so it cannot resolve the topographic
wind acceleration around the Point Conception cape. GWA is a 250 m downscaled
product (mesoscale + microscale modelling over ERA5), fine enough to actually
see the corner-acceleration gradient between e.g. Arguello/Conception (exposed,
windward of the bend) and CDFW (lee, inside the Santa Barbara Channel bend).

Data source: https://globalwindatlas.info (DTU / World Bank / Vortex).
The GIS download API redirects to per-country Cloud-Optimized GeoTIFFs on a
CDN; rather than pulling the whole CONUS raster (tens of GB), we do a
windowed read over just the study-area bounding box via GDAL's /vsicurl/
driver (rasterio), so only the needed tile is transferred.

Only mean wind speed is fetched here (streamlined build). combined-Weibull-A/
-k and power-density are available from the same API at the same resolution
if a wind-speed *distribution* (not just the mean) is wanted later.
"""

import csv
import os

import numpy as np
import rasterio
from rasterio.windows import from_bounds

GWA_URL = "/vsicurl/https://gwa.cdn.nazkamapps.com/country_tifs_v4/USA_wind-speed_10m.tif"

# Study-area bounding box (lon_min, lat_min, lon_max, lat_max), buffered
# ~0.1 deg beyond the site centroids so the clipped raster is still useful
# for a future resistance-surface (rather than just point extraction).
BBOX = (-120.75, 34.35, -120.05, 35.00)

FIELD_DATA = "/private/groups/meyerlab/esteiger/field_data.txt"

OUT_DIR_RAW = "/private/groups/meyerlab/esteiger/wind/gwa"
OUT_TIF = os.path.join(OUT_DIR_RAW, "gwa_windspeed_10m_ptconception.tif")

OUT_DIR_PROJECT = "/private/groups/meyerlab/esteiger/wind"
OUT_SITE_VALUES = os.path.join(OUT_DIR_PROJECT, "gwa_windspeed_site_values.txt")


def compute_centroids(field_data_path):
    """Return, per cluster: centroid (lat, lon) and the sample bounding box
    (lat_min, lat_max, lon_min, lon_max).

    Some clusters (Point_Conception: 14.8 x 7.1 km; Strauss: 2.8 x 2.8 km)
    span an area much larger than one 250 m GWA pixel. A single
    nearest-pixel read at the centroid is noisy in this terrain (Strauss's
    centroid pixel reads 0.56 m/s vs. 4.8-4.9 m/s at Dangermond/Arguello,
    6-10 km away) -- so wind speed is extracted as the mean over each
    cluster's actual sampled footprint, not a single point.
    """
    rows = {}
    with open(field_data_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                lat, lon = float(row["Lat"]), float(row["Long"])
            except ValueError:
                continue
            rows.setdefault(row["Clusters"], []).append((lat, lon))
    centroids = {}
    for cluster, coords in rows.items():
        lats, lons = zip(*coords)
        centroids[cluster] = {
            "lat": sum(lats) / len(lats),
            "lon": sum(lons) / len(lons),
            "lat_min": min(lats), "lat_max": max(lats),
            "lon_min": min(lons), "lon_max": max(lons),
            "n": len(coords),
        }
    return centroids


def main():
    os.makedirs(OUT_DIR_RAW, exist_ok=True)
    os.makedirs(OUT_DIR_PROJECT, exist_ok=True)

    centroids = compute_centroids(FIELD_DATA)
    print(f"Loaded {len(centroids)} population centroids from {FIELD_DATA}")

    print(f"Opening GWA COG (remote, windowed read): {GWA_URL}")
    with rasterio.open(GWA_URL) as src:
        window = from_bounds(*BBOX, transform=src.transform)
        data = src.read(1, window=window)
        clipped_transform = src.window_transform(window)
        profile = src.profile.copy()
        profile.update(
            height=data.shape[0],
            width=data.shape[1],
            transform=clipped_transform,
            driver="GTiff",
            compress="deflate",
        )

        print(f"Clipped raster: {data.shape[1]} x {data.shape[0]} px, "
              f"bbox={BBOX}")
        with rasterio.open(OUT_TIF, "w", **profile) as dst:
            dst.write(data, 1)
        print(f"Saved clipped raster -> {OUT_TIF}")

        # Extract wind speed as the mean over each cluster's sampled
        # footprint (bounding box of its individuals' GPS points, buffered
        # by one pixel ~250 m so small/point-like clusters still get a
        # small neighborhood rather than exactly one raw pixel) from the
        # *original* remote dataset (full precision), not the clipped copy.
        buf = 0.0025  # deg, ~1 GWA pixel
        with open(OUT_SITE_VALUES, "w") as out:
            out.write("Cluster\tLat\tLong\tn_samples\tgwa_windspeed_10m_ms_mean\t"
                      "gwa_windspeed_10m_ms_sd\tn_pixels\n")
            for cluster, c in sorted(centroids.items()):
                win = from_bounds(
                    c["lon_min"] - buf, c["lat_min"] - buf,
                    c["lon_max"] + buf, c["lat_max"] + buf,
                    transform=src.transform,
                )
                vals = src.read(1, window=win).astype(float)
                mean_val, sd_val = float(np.mean(vals)), float(np.std(vals))
                out.write(f"{cluster}\t{c['lat']:.6f}\t{c['lon']:.6f}\t{c['n']}\t"
                          f"{mean_val:.4f}\t{sd_val:.4f}\t{vals.size}\n")
                print(f"  {cluster:20s} n={c['n']:4d}  "
                      f"wind speed = {mean_val:.2f} +/- {sd_val:.2f} m/s "
                      f"(mean over {vals.size} pixels)")

    print(f"Saved site-level values -> {OUT_SITE_VALUES}")


if __name__ == "__main__":
    main()
