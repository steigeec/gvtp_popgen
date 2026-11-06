#!/private/home/esteiger/miniconda3/envs/herbie/bin/python3
"""
Download HRRR 10 m wind (u/v) for the Point Conception study area, bloom
season (Jun-Oct) 2015-2022, 6-hourly (00, 06, 12, 18 UTC).

Input:  HRRR analyses (fxx=0), fetched by Herbie
Output: BASE_DIR/wind/hrrr/hrrr_wind_raw_points.csv
        (datetime, lat, lon, u10, v10). u/v are grid-relative; run
        rotate_hrrr_wind.py next.
Env:    'herbie' conda env
Usage:  /private/home/esteiger/miniconda3/envs/herbie/bin/python3 download_hrrr_wind.py
"""

import os

import pandas as pd
from herbie import FastHerbie

# Buffered study-area bbox; HRRR longitudes are 0-360
LON_MIN, LON_MAX = 239.25, 239.95   # -120.75 to -120.05
LAT_MIN, LAT_MAX = 34.35, 34.95

# 2015-2022: pre-turbine (Strauss Wind Farm began operating 2023)
YEARS = list(range(2015, 2023))
BLOOM_MONTHS = [6, 7, 8, 9, 10]   # Gaviota tarplant flowering season, June-October
HOURS = [0, 6, 12, 18]

BASE_DIR = "/private/groups/meyerlab/esteiger"
OUT_DIR = f"{BASE_DIR}/wind/hrrr"
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


def fetch_datetimes(datetimes):
    """Fetch timesteps in one FastHerbie call. Returns the extracted frames
    and the datetimes that came back missing."""
    try:
        FH = FastHerbie(datetimes, model="hrrr", product="sfc", fxx=[0], max_threads=50)
        ds = FH.xarray(":[U|V]GRD:10 m", max_threads=50)
    except Exception as e:   # one bad file fails the whole chunk, so retry one at a time
        print(f"  chunk failed ({type(e).__name__}: {e})")
        if len(datetimes) == 1:
            return [], datetimes
        frames, missing = [], []
        for dt in datetimes:
            got, miss = fetch_datetimes([dt])
            frames.extend(got)
            missing.extend(miss)
        return frames, missing

    frames, missing = [], []
    for dt in datetimes:
        try:
            sub = ds.sel(time=dt) if "time" in ds.dims else ds
        except KeyError:
            sub = None
        df = extract_bbox(dt, sub)
        if df is not None:
            frames.append(df)
        else:
            missing.append(dt)
    return frames, missing


def fetch_month(year, month, attempts=4):
    """Fetch a month, retrying failed timesteps up to `attempts` times."""
    frames, todo = [], build_datetimes(year, month)
    for _ in range(attempts):
        got, todo = fetch_datetimes(todo)
        frames.extend(got)
        if not todo:
            break
    return frames, len(build_datetimes(year, month)) - len(todo), len(todo)


def main():
    # Resume: months already saved in full are kept and skipped.
    all_frames, done = [], set()
    if os.path.exists(RAW_POINTS_CSV):
        prev = pd.read_csv(RAW_POINTS_CSV, parse_dates=["datetime"])
        ym = list(zip(prev["datetime"].dt.year, prev["datetime"].dt.month))
        n_saved = prev.groupby(pd.Series(ym, index=prev.index))["datetime"].nunique()
        done = {k for k, n in n_saved.items()
                if k[0] in YEARS and k[1] in BLOOM_MONTHS and n == len(build_datetimes(*k))}
        all_frames.append(prev[[k in done for k in ym]])
    total_ok, total_fail = 0, 0
    for year in YEARS:
        for month in BLOOM_MONTHS:
            if (year, month) in done:
                continue
            print(f"=== {year}-{month:02d} ===")
            frames, n_ok, n_fail = fetch_month(year, month)
            all_frames.extend(frames)
            total_ok += n_ok
            total_fail += n_fail
            print(f"  {n_ok} ok, {n_fail} failed/missing "
                  f"(running total: {total_ok} ok, {total_fail} failed)")
            # Written after each month so an interrupted run can resume.
            pd.concat(all_frames, ignore_index=True).to_csv(RAW_POINTS_CSV, index=False)

    print(f"\nSaved raw (u10,v10) points -> {RAW_POINTS_CSV} "
          f"({total_ok} timesteps ok, {total_fail} failed/missing)")


if __name__ == "__main__":
    main()
