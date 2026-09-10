#!/usr/bin/env python3
"""
Parse the raw station wind observations downloaded by download_station_wind.sh
into one tidy table (station, datetime, dir_deg, speed_ms), then summarize
each station as a 16-sector wind rose (direction-bin frequency + mean speed)
and a circular resultant vector (prevailing direction + steadiness).

NDBC's historical file format changed over the 1995-2023 PTGC1 record
(column set/header rows differ by era), so columns are located by name from
each file's own header rather than assumed to be in a fixed position.

Calm/missing handling:
  - NDBC: WD/WDIR == 999 is the "missing or variable direction" sentinel.
  - ASOS (IEM): blank drct is missing; dir == 0 & speed == 0 is calm
    (by convention -- true north is coded 360, not 0).
Both are excluded from the direction-binned wind rose but calm fraction is
reported, since a high calm fraction matters for interpreting wind-driven
dispersal potential.
"""

import glob
import gzip
import math
import os

import numpy as np
import pandas as pd

STATION_OBS_DIR = "/private/groups/meyerlab/esteiger/wind/station_obs"
OUT_DIR = "/private/groups/meyerlab/esteiger/wind"
# Combined tidy observation table is raw-ish (no analysis applied) and large
# (~15 MB) -- lives with the other raw downloads in env/wind/, not in the
# git-tracked project directory. windrose_summary/sectors ARE derived
# analysis outputs (small) and stay in the project dir with the scripts.
OUT_OBS = "/private/groups/meyerlab/esteiger/wind/station_wind_obs.csv"
OUT_ROSE = os.path.join(OUT_DIR, "windrose_summary.txt")
OUT_SECTORS = os.path.join(OUT_DIR, "windrose_sectors.txt")

KT_TO_MS = 0.514444
N_SECTORS = 16
SECTOR_WIDTH = 360 / N_SECTORS

# Wind is hypothesized to matter here via pollen/seed transport, which can
# only happen while the plant is actually flowering/dispersing -- June
# through September. The wind-rose statistics (prevailing direction, speed,
# steadiness) are therefore computed on bloom-season observations only, not
# an annual average, since off-season wind (e.g. winter storm patterns) is
# not mechanistically relevant to the hypothesis being tested.
BLOOM_MONTHS = [6, 7, 8, 9]


def parse_ptgc1():
    rows = []
    files = sorted(glob.glob(os.path.join(STATION_OBS_DIR, "ptgc1", "*.txt.gz")))
    for f in files:
        with gzip.open(f, "rt") as fh:
            lines = fh.readlines()
        header = lines[0].lstrip("#").split()
        data_start = 1
        if lines[1].startswith("#"):
            data_start = 2
        dir_col = "WDIR" if "WDIR" in header else "WD"
        dir_idx = header.index(dir_col)
        spd_idx = header.index("WSPD")
        has_mm = "mm" in header or "MM" in [h.upper() for h in header]
        yr_idx, mo_idx, dy_idx, hh_idx = 0, 1, 2, 3
        mn_idx = 4 if header[4].lower() == "mm" else None
        for line in lines[data_start:]:
            parts = line.split()
            if len(parts) <= max(dir_idx, spd_idx):
                continue
            try:
                yr = int(parts[yr_idx])
                yr = yr + 1900 if yr < 100 else yr
                mo, dy, hh = int(parts[mo_idx]), int(parts[dy_idx]), int(parts[hh_idx])
                mn = int(parts[mn_idx]) if mn_idx is not None else 0
                dir_deg = float(parts[dir_idx])
                speed = float(parts[spd_idx])  # already m/s in NDBC files
            except (ValueError, IndexError):
                continue
            if dir_deg == 999:
                dir_deg = np.nan
            rows.append(("PTGC1", f"{yr:04d}-{mo:02d}-{dy:02d} {hh:02d}:{mn:02d}",
                         dir_deg, speed))
    return pd.DataFrame(rows, columns=["station", "datetime", "dir_deg", "speed_ms"])


def parse_asos(station):
    f = os.path.join(STATION_OBS_DIR, f"{station}_asos_2000_2024.csv")
    df = pd.read_csv(f, na_values=[""])
    df = df.rename(columns={"valid": "datetime", "drct": "dir_deg", "sknt": "speed_ms"})
    df["speed_ms"] = df["speed_ms"] * KT_TO_MS
    # calm convention: dir==0 & speed==0 -> calm, not "due north"
    calm = (df["dir_deg"] == 0) & (df["speed_ms"] == 0)
    df.loc[calm, "dir_deg"] = np.nan
    return df[["station", "datetime", "dir_deg", "speed_ms"]]


def circular_resultant(dir_deg, speed):
    """Mean wind vector (speed-weighted): resultant direction + resultant
    length R in [0,1] (1 = perfectly steady direction, 0 = fully variable)."""
    rad = np.radians(dir_deg.values)
    u = -speed.values * np.sin(rad)  # eastward component (met convention: dir = FROM)
    v = -speed.values * np.cos(rad)  # northward component
    u_mean, v_mean = np.nanmean(u), np.nanmean(v)
    mean_speed = np.nanmean(speed.values)
    resultant_speed = math.hypot(u_mean, v_mean)
    resultant_dir = (math.degrees(math.atan2(-u_mean, -v_mean))) % 360
    R = resultant_speed / mean_speed if mean_speed > 0 else np.nan
    return resultant_dir, R, mean_speed


def sector_table(sub):
    valid_dir = sub.dropna(subset=["dir_deg"])
    n_total = len(sub)
    n_calm_missing = n_total - len(valid_dir)
    sector = (np.floor((valid_dir["dir_deg"] + SECTOR_WIDTH / 2) / SECTOR_WIDTH)
              .astype(int) % N_SECTORS)
    out = []
    for s in range(N_SECTORS):
        center = s * SECTOR_WIDTH
        in_sector = valid_dir[sector == s]
        freq = len(in_sector) / n_total if n_total else np.nan
        mean_spd = in_sector["speed_ms"].mean() if len(in_sector) else 0.0
        out.append((center, freq, mean_spd, len(in_sector)))
    df = pd.DataFrame(out, columns=["sector_center_deg", "frequency", "mean_speed_ms", "n_obs"])
    df["calm_or_missing_frac"] = n_calm_missing / n_total if n_total else np.nan
    return df


def main():
    ptgc1 = parse_ptgc1()
    print(f"PTGC1: {len(ptgc1)} raw observations parsed")

    frames = [ptgc1]
    for st in ["VGN", "TQS", "SBA"]:
        df = parse_asos(st)
        print(f"{st}: {len(df)} raw observations parsed")
        frames.append(df)

    obs = pd.concat(frames, ignore_index=True)
    obs.to_csv(OUT_OBS, index=False)
    print(f"Saved combined observation table (all months) -> {OUT_OBS} ({len(obs)} rows)")

    obs_month = pd.to_datetime(obs["datetime"]).dt.month
    bloom_obs = obs[obs_month.isin(BLOOM_MONTHS)].copy()
    print(f"Bloom-season (Jun-Sep) subset: {len(bloom_obs)} / {len(obs)} observations "
          f"({100 * len(bloom_obs) / len(obs):.1f}%)")
    print("Wind-rose statistics below are computed on the bloom-season subset only.")

    rose_rows = []
    sector_frames = []
    for station, sub in bloom_obs.groupby("station"):
        n_total = len(sub)
        n_valid_dir = sub["dir_deg"].notna().sum()
        resultant_dir, R, mean_speed = circular_resultant(
            sub.dropna(subset=["dir_deg", "speed_ms"])["dir_deg"],
            sub.dropna(subset=["dir_deg", "speed_ms"])["speed_ms"],
        )
        rose_rows.append((station, n_total, n_valid_dir,
                           round(mean_speed, 3), round(resultant_dir, 1), round(R, 3)))
        sec = sector_table(sub)
        sec.insert(0, "station", station)
        sector_frames.append(sec)

    rose = pd.DataFrame(rose_rows, columns=[
        "station", "n_obs", "n_valid_dir", "mean_speed_ms",
        "resultant_dir_deg_from", "resultant_R"])
    rose.to_csv(OUT_ROSE, sep="\t", index=False)
    print(f"\nSaved wind-rose summary -> {OUT_ROSE}")
    print(rose.to_string(index=False))

    sectors = pd.concat(sector_frames, ignore_index=True)
    sectors.to_csv(OUT_SECTORS, sep="\t", index=False)
    print(f"Saved 16-sector tables -> {OUT_SECTORS}")


if __name__ == "__main__":
    main()
