#!/private/home/esteiger/miniconda3/envs/herbie/bin/python3
"""
Rotate the raw HRRR 10 m u/v (download_hrrr_wind.py) from grid-relative to
earth-relative (true east/north).

HRRR GRIB2 u/v are relative to the Lambert Conformal grid axes
(uvRelativeToGrid = 1), and neither Herbie nor cfgrib rotates them. The
rotation angle is

    alpha = sin(Latin1) * (lon - LoV)

with Latin1 = 38.5 and LoV = 262.5 deg (about -14.3 deg at the study area), and

    u_earth =  cos(alpha) * u_grid + sin(alpha) * v_grid
    v_earth = -sin(alpha) * u_grid + cos(alpha) * v_grid

Cell coordinates differ by ~0.0006 deg between HRRR versions, so all points
are snapped to the latest timestep's coordinates.

Input:  BASE_DIR/wind/hrrr/hrrr_wind_raw_points.csv
Output: BASE_DIR/wind/hrrr/hrrr_wind_raw_points_earth.csv
        (same columns; used by all downstream scripts)
Env:    'herbie' conda env
Usage:  /private/home/esteiger/miniconda3/envs/herbie/bin/python3 rotate_hrrr_wind.py
"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

BASE_DIR = "/private/groups/meyerlab/esteiger"
HRRR_DIR = f"{BASE_DIR}/wind/hrrr"
IN_CSV = f"{HRRR_DIR}/hrrr_wind_raw_points.csv"
OUT_CSV = f"{HRRR_DIR}/hrrr_wind_raw_points_earth.csv"

LATIN1 = 38.5   # HRRR Lambert standard parallel (deg)
LOV = 262.5     # HRRR Lambert central meridian (deg, 0-360)


def main():
    df = pd.read_csv(IN_CSV)
    ref = df.loc[df["datetime"] == df["datetime"].max(), ["lat", "lon"]].values
    dist, idx = cKDTree(ref).query(df[["lat", "lon"]].values)
    assert dist.max() < 0.005, f"cell coordinate offset {dist.max():.4f} deg"
    df["lat"], df["lon"] = ref[idx, 0], ref[idx, 1]
    assert not df.duplicated(["datetime", "lat", "lon"]).any()
    dlon = (df["lon"] - LOV + 180) % 360 - 180   # lon is 0-360, like LoV
    alpha = np.radians(np.sin(np.radians(LATIN1)) * dlon)
    u_g, v_g = df["u10"].values, df["v10"].values
    df["u10"] = np.cos(alpha) * u_g + np.sin(alpha) * v_g
    df["v10"] = -np.sin(alpha) * u_g + np.cos(alpha) * v_g

    print(f"Rotation angle across {df[['lat', 'lon']].drop_duplicates().shape[0]} grid cells: "
          f"{np.degrees(alpha).min():.2f} to {np.degrees(alpha).max():.2f} deg")
    speed_change = np.abs(np.hypot(df["u10"], df["v10"]) - np.hypot(u_g, v_g)).max()
    print(f"Max |speed change| (should be ~0): {speed_change:.2e} m/s")

    df.to_csv(OUT_CSV, index=False)
    print(f"Saved earth-relative u/v -> {OUT_CSV} ({len(df)} rows)")


if __name__ == "__main__":
    main()
