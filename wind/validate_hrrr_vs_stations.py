#!/private/home/esteiger/miniconda3/envs/herbie/bin/python3
"""
Compare raw HRRR wind (earth-relative) at the nearest grid cell to each
station against the station observations, over all bloom-season timesteps
from download_hrrr_wind.py.

Stations (coordinates as in download_station_wind.sh):
  PTGC1 (34.577, -120.648)   NDBC C-MAN, Point Arguello
  VGN   (34.8495, -120.5996) Vandenberg SFB North ASOS
  TQS   (34.636, -120.6161)  Vandenberg SFB South ASOS
Only PTGC1 has observations before 2023.

Input:  BASE_DIR/wind/hrrr/hrrr_wind_raw_points_earth.csv,
        BASE_DIR/wind/station_wind_obs.csv (PTGC1),
        BASE_DIR/wind/station_obs/{VGN,TQS}_asos_2000_2024.csv
Output: BASE_DIR/wind/hrrr/validate_{station}_matched.csv,
        BASE_DIR/wind/hrrr_vs_station_validation_summary.csv
Env:    'herbie' conda env
Usage:  /private/home/esteiger/miniconda3/envs/herbie/bin/python3 validate_hrrr_vs_stations.py
"""

import numpy as np
import pandas as pd

BASE_DIR = "/private/groups/meyerlab/esteiger"
WIND_DIR = f"{BASE_DIR}/wind"
HRRR_CSV = f"{WIND_DIR}/hrrr/hrrr_wind_raw_points_earth.csv"
PTGC1_CSV = f"{WIND_DIR}/station_wind_obs.csv"

STATIONS = {
    "PTGC1": dict(lat=34.577, lon=-120.648, kind="buoy"),
    "VGN":   dict(lat=34.8495, lon=-120.5996, kind="asos"),
    "TQS":   dict(lat=34.636, lon=-120.6161, kind="asos"),
}

TIME_TOL = pd.Timedelta("35min")  # ASOS obs land near :55, not on the hour


def to_from_direction(u, v):
    """Meteorological direction the wind blows from, in compass degrees."""
    return (np.degrees(np.arctan2(-u, -v))) % 360


def circular_abs_error(a, b):
    d = np.abs(a - b) % 360
    return np.minimum(d, 360 - d)


def load_hrrr():
    df = pd.read_csv(HRRR_CSV, parse_dates=["datetime"])
    df["lon180"] = df["lon"] - 360
    return df


def nearest_cell_series(hrrr, lat, lon):
    """Time series of the HRRR cell nearest to (lat, lon)."""
    first_t = hrrr["datetime"].iloc[0]
    snap = hrrr[hrrr["datetime"] == first_t]
    d = np.hypot(snap["lat"] - lat, snap["lon180"] - lon)
    cell_lat, cell_lon = snap.loc[d.idxmin(), ["lat", "lon180"]]
    dist_km = d.min() * 111.0  # rough degrees to km at this latitude
    sub = hrrr[(hrrr["lat"] == cell_lat) & (hrrr["lon180"] == cell_lon)].copy()
    sub["speed_hrrr"] = np.hypot(sub["u10"], sub["v10"])
    sub["dir_hrrr"] = to_from_direction(sub["u10"], sub["v10"])
    return sub[["datetime", "speed_hrrr", "dir_hrrr"]], dist_km


def load_ptgc1_obs():
    df = pd.read_csv(PTGC1_CSV, parse_dates=["datetime"])
    df = df[df["station"] == "PTGC1"].copy()
    df = df.rename(columns={"dir_deg": "dir_obs", "speed_ms": "speed_obs"})
    return df[["datetime", "dir_obs", "speed_obs"]]


def load_asos_obs(station):
    df = pd.read_csv(f"{WIND_DIR}/station_obs/{station}_asos_2000_2024.csv")
    df["datetime"] = pd.to_datetime(df["valid"])
    df = df.rename(columns={"drct": "dir_obs"})
    df["speed_obs"] = pd.to_numeric(df["sknt"], errors="coerce") * 0.514444
    df["dir_obs"] = pd.to_numeric(df["dir_obs"], errors="coerce")
    return df[["datetime", "dir_obs", "speed_obs"]].dropna(subset=["speed_obs"])


def match_nearest_time(hrrr_sub, obs, tol=TIME_TOL):
    hrrr_sub = hrrr_sub.sort_values("datetime")
    obs = obs.sort_values("datetime")
    merged = pd.merge_asof(
        hrrr_sub, obs, on="datetime", direction="nearest", tolerance=tol
    )
    return merged.dropna(subset=["speed_obs"])


def summarize(name, dist_km, m):
    n = len(m)
    if n == 0:
        print(f"\n=== {name} (nearest HRRR cell {dist_km:.1f} km away) ===")
        print("  no matching timesteps")
        return
    speed_bias = (m["speed_obs"] - m["speed_hrrr"]).mean()
    speed_mae = (m["speed_obs"] - m["speed_hrrr"]).abs().mean()
    speed_rmse = np.sqrt(((m["speed_obs"] - m["speed_hrrr"]) ** 2).mean())
    speed_corr = m["speed_obs"].corr(m["speed_hrrr"])

    calm = (m["speed_obs"] < 0.5) | (m["speed_hrrr"] < 0.5)
    m_dir = m[~calm]
    dir_mae = circular_abs_error(m_dir["dir_obs"], m_dir["dir_hrrr"]).mean()
    dir_median = circular_abs_error(m_dir["dir_obs"], m_dir["dir_hrrr"]).median()

    # vector correlation: mean of the u and v correlations
    u_obs = -m_dir["speed_obs"] * np.sin(np.radians(m_dir["dir_obs"]))
    v_obs = -m_dir["speed_obs"] * np.cos(np.radians(m_dir["dir_obs"]))
    u_hrrr = -m_dir["speed_hrrr"] * np.sin(np.radians(m_dir["dir_hrrr"]))
    v_hrrr = -m_dir["speed_hrrr"] * np.cos(np.radians(m_dir["dir_hrrr"]))
    vec_corr = (u_obs.corr(u_hrrr) + v_obs.corr(v_hrrr)) / 2

    years = sorted(m["datetime"].dt.year.unique())
    print(f"\n=== {name} (nearest HRRR cell {dist_km:.1f} km away) ===")
    print(f"  n = {n} matched timesteps, years covered: {years}")
    print(f"  speed: obs mean {m['speed_obs'].mean():.2f} m/s, "
          f"HRRR mean {m['speed_hrrr'].mean():.2f} m/s, "
          f"bias(obs-HRRR) {speed_bias:+.2f}, MAE {speed_mae:.2f}, "
          f"RMSE {speed_rmse:.2f}, corr {speed_corr:.2f}")
    print(f"  direction (excl. {calm.sum()} calm obs): "
          f"MAE {dir_mae:.0f} deg, median {dir_median:.0f} deg, "
          f"u/v vector corr {vec_corr:.2f}")
    return dict(station=name, n=n, dist_km=dist_km, speed_bias=speed_bias,
                speed_mae=speed_mae, speed_rmse=speed_rmse,
                speed_corr=speed_corr, dir_mae=dir_mae, vec_corr=vec_corr)


def main():
    hrrr = load_hrrr()
    results = []

    for name, meta in STATIONS.items():
        hrrr_sub, dist_km = nearest_cell_series(hrrr, meta["lat"], meta["lon"])
        obs = load_ptgc1_obs() if name == "PTGC1" else load_asos_obs(name)
        m = match_nearest_time(hrrr_sub, obs)
        r = summarize(name, dist_km, m)
        if r:
            results.append(r)
        out_csv = f"{WIND_DIR}/hrrr/validate_{name}_matched.csv"
        m.to_csv(out_csv, index=False)
        print(f"  matched pairs -> {out_csv}")

    if results:
        summary = pd.DataFrame(results)
        summary_csv = f"{WIND_DIR}/hrrr_vs_station_validation_summary.csv"
        summary.to_csv(summary_csv, index=False)
        print(f"\nSummary table -> {summary_csv}")
        print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
