#!/private/home/esteiger/miniconda3/envs/herbie/bin/python3
"""
Compare WindNinja-downscaled wind at PTGC1's own coordinates against raw
HRRR and PTGC1 observations, across the 486 matched 2022 bloom-season
timesteps -- the direct test of whether terrain downscaling closes the
~28% raw-HRRR speed underestimate found at this exposed cape station
(see hrrr_vs_station_validation_summary.csv).
"""
import glob
import re
import subprocess

import numpy as np
import pandas as pd

WORKDIR = "/private/groups/meyerlab/esteiger/wind/windninja_ptgc1"
PTGC1_LAT, PTGC1_LON = 34.577, -120.648

BASE = pd.read_csv(
    "/private/groups/meyerlab/esteiger/wind/hrrr/validate_PTGC1_matched.csv",
    parse_dates=["datetime"],
)


def sample(asc_path, lat, lon):
    out = subprocess.run(
        ["gdallocationinfo", "-wgs84", "-valonly", asc_path, str(lon), str(lat)],
        capture_output=True, text=True,
    )
    v = out.stdout.strip()
    return float(v) if v else np.nan


def parse_dt_from_name(path):
    # e.g. ptgc1_dem_08-26-2022_1700_460m_vel.asc -> local time embedded
    m = re.search(r"_(\d{2})-(\d{2})-(\d{4})_(\d{2})(\d{2})_460m_vel\.asc$", path)
    mo, d, y, hh, mm = m.groups()
    local = pd.Timestamp(int(y), int(mo), int(d), int(hh), int(mm), tz="America/Los_Angeles")
    return local.tz_convert("UTC").tz_localize(None)


def circular_abs_error(a, b):
    d = np.abs(a - b) % 360
    return np.minimum(d, 360 - d)


def main():
    vel_files = glob.glob(f"{WORKDIR}/grib2/hrrr/*/*_vel.asc")
    print(f"Found {len(vel_files)} WindNinja output files (expected 486)")

    rows = []
    for vf in vel_files:
        af = vf.replace("_vel.asc", "_ang.asc")
        dt = parse_dt_from_name(vf)
        speed = sample(vf, PTGC1_LAT, PTGC1_LON)
        direction = sample(af, PTGC1_LAT, PTGC1_LON)
        rows.append(dict(datetime=dt, speed_wn=speed, dir_wn=direction))

    wn = pd.DataFrame(rows)
    m = BASE.merge(wn, on="datetime", how="inner")
    print(f"Merged {len(m)} of {len(BASE)} matched timesteps")

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

        print(f"\n=== {label} vs PTGC1 obs (n={len(m)}) ===")
        print(f"  obs mean {m['speed_obs'].mean():.2f} m/s, {label} mean {m[speed_col].mean():.2f} m/s")
        print(f"  bias(obs-{label}) {bias:+.2f} m/s ({pct_bias:+.0f}% of obs mean), "
              f"MAE {mae:.2f}, RMSE {rmse:.2f}, corr {corr:.2f}")
        print(f"  direction MAE {dir_mae:.0f} deg, median {dir_median:.0f} deg")

    out_csv = f"{WORKDIR}/ptgc1_hrrr_vs_windninja_vs_obs.csv"
    m.to_csv(out_csv, index=False)
    print(f"\nFull comparison table -> {out_csv}")


if __name__ == "__main__":
    main()
