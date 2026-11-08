#!/private/home/esteiger/miniconda3/envs/herbie/bin/python3
"""
Compare WindNinja wind at each station's grid cell (PTGC1, VGN, TQS) against
raw HRRR and the station observations, for the timesteps already matched by
validate_hrrr_vs_stations.py. Uses the production run from
run_windninja_landfire_batch.sh. Only PTGC1 has observations before 2023, so
VGN and TQS have no matched timesteps and are skipped.

Input:  BASE_DIR/wind/windninja_landfire/region/*_vel.asc, *_ang.asc
        BASE_DIR/wind/hrrr/validate_{PTGC1,VGN,TQS}_matched.csv
Output: windninja_landfire/region/{station}_hrrr_vs_windninja_vs_obs.csv,
        BASE_DIR/wind/windninja_station_validation_summary.csv
Conda env: herbie
  conda create -n herbie -c conda-forge python=3.11 herbie-data=2026.3.0
  (Herbie also brings pandas, xarray, cfgrib and scipy)
Also needs gdallocationinfo on PATH (from the windninja env's GDAL).
Usage:  PATH=/private/home/esteiger/miniconda3/envs/windninja/bin:$PATH \
          /private/home/esteiger/miniconda3/envs/herbie/bin/python3 validate_windninja_stations.py
"""
import glob
import re
import subprocess

import numpy as np
import pandas as pd

BASE_DIR = "/private/groups/meyerlab/esteiger"
WIND_DIR = f"{BASE_DIR}/wind"
WORKDIR = f"{WIND_DIR}/windninja_landfire/region"

STATIONS = {
    "PTGC1": dict(lat=34.577, lon=-120.648),
    "VGN": dict(lat=34.8495, lon=-120.5996),
    "TQS": dict(lat=34.636, lon=-120.6161),
}


def station_cell(asc_path, lat, lon):
    """(row, col) of the grid cell containing the station; all outputs share one grid."""
    out = subprocess.run(["gdallocationinfo", "-wgs84", "-xml", asc_path, str(lon), str(lat)],
                         capture_output=True, text=True, check=True).stdout
    m = re.search(r'pixel="(\d+)" line="(\d+)"', out)
    return int(m.group(2)), int(m.group(1))


def sample(asc_path, cell):
    v = np.loadtxt(asc_path, skiprows=6)[cell]
    return np.nan if v == -9999 else float(v)


def parse_dt_from_name(path):
    # WindNinja output name, e.g. region_landfire_05-31-2022_1700_460m_non_neutral_stability_vel.asc
    # The embedded time is local (America/Los_Angeles).
    m = re.search(r"_(\d{2})-(\d{2})-(\d{4})_(\d{2})(\d{2})_460m_non_neutral_stability_vel\.asc$", path)
    mo, d, y, hh, mm = m.groups()
    local = pd.Timestamp(int(y), int(mo), int(d), int(hh), int(mm), tz="America/Los_Angeles")
    return local.tz_convert("UTC").tz_localize(None)


def circular_abs_error(a, b):
    d = np.abs(a - b) % 360
    return np.minimum(d, 360 - d)


def main():
    vel_files = glob.glob(f"{WORKDIR}/*_vel.asc")
    print(f"Found {len(vel_files)} WindNinja output files (4741 HRRR timesteps)")
    vel_by_dt = {parse_dt_from_name(vf): vf for vf in vel_files}

    summary_rows = []
    for name, meta in STATIONS.items():
        base = pd.read_csv(f"{WIND_DIR}/hrrr/validate_{name}_matched.csv", parse_dates=["datetime"])
        if base.empty:
            print(f"\n##### {name}: no observations in the wind years, skipped")
            continue

        cell = station_cell(vel_files[0], meta["lat"], meta["lon"])
        rows = []
        for dt in base["datetime"]:
            vf = vel_by_dt.get(dt)
            if vf is None:
                continue
            af = vf.replace("_vel.asc", "_ang.asc")
            rows.append(dict(datetime=dt,
                             speed_wn=sample(vf, cell),
                             dir_wn=sample(af, cell)))
        wn = pd.DataFrame(rows)
        m = base.merge(wn, on="datetime", how="inner").dropna(subset=["speed_wn", "dir_wn"])
        print(f"\n##### {name}: merged {len(m)} of {len(base)} matched timesteps")

        for label, speed_col, dir_col in [
            ("raw HRRR", "speed_hrrr", "dir_hrrr"),
            ("WindNinja", "speed_wn", "dir_wn"),
        ]:
            bias = (m["speed_obs"] - m[speed_col]).mean()
            mae = (m["speed_obs"] - m[speed_col]).abs().mean()
            rmse = np.sqrt(((m["speed_obs"] - m[speed_col]) ** 2).mean())
            corr = m["speed_obs"].corr(m[speed_col])
            pct_bias = 100 * bias / m["speed_obs"].mean()

            calm = (m["speed_obs"] < 0.5) | (m[speed_col] < 0.5)
            dsub = m[~calm]
            dir_mae = circular_abs_error(dsub["dir_obs"], dsub[dir_col]).mean()
            dir_median = circular_abs_error(dsub["dir_obs"], dsub[dir_col]).median()

            print(f"\n=== {label} vs {name} obs (n={len(m)}) ===")
            print(f"  obs mean {m['speed_obs'].mean():.2f} m/s, {label} mean {m[speed_col].mean():.2f} m/s")
            print(f"  bias(obs-{label}) {bias:+.2f} m/s ({pct_bias:+.0f}% of obs mean), "
                  f"MAE {mae:.2f}, RMSE {rmse:.2f}, corr {corr:.2f}")
            print(f"  direction MAE {dir_mae:.0f} deg, median {dir_median:.0f} deg")
            summary_rows.append(dict(station=name, product=label, n=len(m),
                                     speed_bias=bias, speed_pct_bias=pct_bias,
                                     speed_mae=mae, speed_rmse=rmse, speed_corr=corr,
                                     n_dir=len(dsub), dir_mae=dir_mae, dir_median=dir_median))

        out_csv = f"{WORKDIR}/{name.lower()}_hrrr_vs_windninja_vs_obs.csv"
        m.to_csv(out_csv, index=False)
        print(f"\nFull comparison table -> {out_csv}")

    summary = pd.DataFrame(summary_rows)
    summary_csv = f"{WIND_DIR}/windninja_station_validation_summary.csv"
    summary.to_csv(summary_csv, index=False)
    print(f"\nSummary table -> {summary_csv}")
    print(summary.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
