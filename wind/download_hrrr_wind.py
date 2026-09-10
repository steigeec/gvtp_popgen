#!/private/home/esteiger/miniconda3/envs/herbie/bin/python3
"""
Download HRRR 10 m wind (U/V components) for the Point Conception study area
and build a bloom-season (Jun-Sep) wind climatology that has BOTH speed and
direction as a genuinely gridded field -- unlike the Global Wind Atlas layer
(speed only, no direction) and the 4 weather stations (direction only at
those 4 points; see chat discussion).

HRRR is a 3 km CONUS model, hourly, back to ~2014, on a curvilinear (Lambert
Conformal) grid -- NOT a simple lat/lon grid, so subsetting is done by
boolean masking on the 2D lat/lon coordinate arrays rather than simple index
slicing.

Sampling scheme: 3 recent years (2022-2024) x bloom months (Jun-Sep) x
6-hourly synoptic times (00, 06, 12, 18 UTC) = ~1460 timesteps. This is a
climatology sample, not exhaustive hourly coverage -- a reasonable trade-off
between capturing the diurnal cycle (6-hourly, standard synoptic sampling)
and download/runtime cost. Extend YEARS / add more hours if a finer
climatology is wanted later.

Run with the isolated 'herbie' conda env (not the shared system/user python):
  /private/home/esteiger/miniconda3/envs/herbie/bin/python3 download_hrrr_wind.py
"""

import os

import pandas as pd
from herbie import FastHerbie

# Buffered bbox matching the GWA download (0-360 longitude convention for HRRR)
LON_MIN, LON_MAX = 239.25, 239.95   # -120.75 to -120.05
LAT_MIN, LAT_MAX = 34.35, 34.95

YEARS = [2022, 2023, 2024]
BLOOM_MONTHS = [6, 7, 8, 9]
HOURS = [0, 6, 12, 18]

OUT_DIR = "/private/groups/meyerlab/esteiger/wind/hrrr"
os.makedirs(OUT_DIR, exist_ok=True)
RAW_POINTS_CSV = os.path.join(OUT_DIR, "hrrr_wind_raw_points.csv")


def build_datetimes(year, month):
    days_in_month = pd.Period(f"{year}-{month:02d}").days_in_month
    return [pd.Timestamp(year=year, month=month, day=day, hour=hour)
            for day in range(1, days_in_month + 1) for hour in HOURS]


def extract_bbox(dt, ds):
    if ds is None:
        return None
    lat = ds["latitude"].values
    lon = ds["longitude"].values
    mask = (lon >= LON_MIN) & (lon <= LON_MAX) & (lat >= LAT_MIN) & (lat <= LAT_MAX)
    if not mask.any():
        return None
    return pd.DataFrame({
        "datetime": dt,
        "lat": lat[mask], "lon": lon[mask],
        "u10": ds["u10"].values[mask], "v10": ds["v10"].values[mask],
    })


def fetch_month(year, month):
    """One FastHerbie call per (year, month) chunk (~120 timesteps) -- keeps
    each call's combined full-CONUS xarray Dataset (~2 GB/chunk) bounded,
    rather than materializing all ~1460 timesteps (~25+ GB) in one merge."""
    datetimes = build_datetimes(year, month)
    FH = FastHerbie(datetimes, model="hrrr", product="sfc", fxx=[0], max_threads=50)
    ds = FH.xarray(":[U|V]GRD:10 m", max_threads=50)

    frames = []
    n_ok, n_fail = 0, 0
    for dt in datetimes:
        try:
            sub = ds.sel(time=dt) if "time" in ds.dims else ds
        except KeyError:
            sub = None
        df = extract_bbox(dt, sub)
        if df is not None:
            frames.append(df)
            n_ok += 1
        else:
            n_fail += 1
    return frames, n_ok, n_fail


def main():
    all_frames = []
    total_ok, total_fail = 0, 0
    for year in YEARS:
        for month in BLOOM_MONTHS:
            print(f"=== {year}-{month:02d} ===")
            frames, n_ok, n_fail = fetch_month(year, month)
            all_frames.extend(frames)
            total_ok += n_ok
            total_fail += n_fail
            print(f"  {n_ok} ok, {n_fail} failed/missing "
                  f"(running total: {total_ok} ok, {total_fail} failed)")
            # Write incrementally so a partial run still leaves usable output.
            pd.concat(all_frames, ignore_index=True).to_csv(RAW_POINTS_CSV, index=False)

    print(f"\nSaved raw (u10,v10) points -> {RAW_POINTS_CSV} "
          f"({total_ok} timesteps ok, {total_fail} failed/missing)")


if __name__ == "__main__":
    main()
