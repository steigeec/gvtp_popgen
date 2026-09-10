#!/usr/bin/env python3
"""
Turn the raw HRRR (u10, v10) points (download_hrrr_wind.py) into a
bloom-season wind climatology, spatially continuous across the whole study
area (478 HRRR 3 km grid cells) rather than limited to the 4 weather
stations: per-grid-cell mean speed, resultant ("prevailing") direction,
directional steadiness (R), and a 16-sector wind-rose table per cell -- the
same quantities build_windrose_table.py computes per station, just computed
per HRRR grid cell instead so every population can use its own local wind
rose rather than borrowing one from a distant station.

Sign convention note: HRRR u10/v10 are true velocity components (TOWARD
convention: +u = blowing toward east, +v = blowing toward north) -- the
opposite convention from station dir_deg, which reports the compass bearing
the wind blows FROM. Everything here is converted to "from" bearings so it's
directly comparable with the station-based windrose_sectors.txt /
windrose_summary.txt already used elsewhere (station_map in
wind_connectivity.R, wind_support()'s half-cosine kernel).
"""

import os

import numpy as np
import pandas as pd

RAW_CSV = "/private/groups/meyerlab/esteiger/wind/hrrr/hrrr_wind_raw_points.csv"
OUT_DIR = "/private/groups/meyerlab/esteiger/wind/hrrr"
OUT_SUMMARY = os.path.join(OUT_DIR, "hrrr_windrose_summary.txt")
OUT_SECTORS = os.path.join(OUT_DIR, "hrrr_windrose_sectors.txt")

N_SECTORS = 16
SECTOR_WIDTH = 360 / N_SECTORS


def to_from_bearing(u, v):
    """HRRR u,v (TOWARD convention) -> compass bearing wind blows FROM."""
    dir_to = np.degrees(np.arctan2(u, v)) % 360
    return (dir_to + 180) % 360


def main():
    df = pd.read_csv(RAW_CSV)
    df["speed"] = np.hypot(df["u10"], df["v10"])
    df["dir_from"] = to_from_bearing(df["u10"].values, df["v10"].values)

    # Grid cell = unique (lat, lon) pair (fixed HRRR grid, exact recurrence
    # across timesteps -- 478 cells x 1464 timesteps = 699792 rows, checked).
    df["cell_id"] = df["lat"].round(6).astype(str) + "_" + df["lon"].round(6).astype(str)

    summary_rows = []
    sector_frames = []
    for cell_id, sub in df.groupby("cell_id"):
        lat0, lon0 = sub["lat"].iloc[0], sub["lon"].iloc[0]
        mean_u, mean_v = sub["u10"].mean(), sub["v10"].mean()
        mean_speed = sub["speed"].mean()
        resultant_speed = np.hypot(mean_u, mean_v)
        R = resultant_speed / mean_speed if mean_speed > 0 else np.nan
        resultant_dir = to_from_bearing(np.array([mean_u]), np.array([mean_v]))[0]

        summary_rows.append((cell_id, lat0, lon0, len(sub),
                              round(mean_speed, 3), round(resultant_dir, 1), round(R, 3)))

        sector = np.floor((sub["dir_from"] + SECTOR_WIDTH / 2) / SECTOR_WIDTH).astype(int) % N_SECTORS
        for s in range(N_SECTORS):
            in_sector = sub[sector == s]
            freq = len(in_sector) / len(sub)
            mean_spd_sector = in_sector["speed"].mean() if len(in_sector) else 0.0
            sector_frames.append((cell_id, s * SECTOR_WIDTH, freq, mean_spd_sector, len(in_sector)))

    summary = pd.DataFrame(summary_rows, columns=[
        "cell_id", "lat", "lon", "n_obs", "mean_speed_ms",
        "resultant_dir_deg_from", "resultant_R"])
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"Saved per-gridcell wind-rose summary -> {OUT_SUMMARY} ({len(summary)} cells)")
    print(summary.describe().to_string())

    sectors = pd.DataFrame(sector_frames, columns=[
        "cell_id", "sector_center_deg", "frequency", "mean_speed_ms", "n_obs"])
    sectors.to_csv(OUT_SECTORS, sep="\t", index=False)
    print(f"Saved per-gridcell 16-sector tables -> {OUT_SECTORS} ({len(sectors)} rows)")


if __name__ == "__main__":
    main()
