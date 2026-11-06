#!/usr/bin/env python3
"""
Parse the raw station wind observations (download_station_wind.sh) into one
table with columns station, datetime, dir_deg, speed_ms.

NDBC column layout changes between years, so PTGC1 columns are found by name
from each file's header. Missing or calm direction is stored as NaN: NDBC
uses 999, and for ASOS dir 0 with speed 0 means calm (north is coded 360).

Input:  BASE_DIR/wind/station_obs/
Output: BASE_DIR/wind/station_wind_obs.csv
Usage:  python3 build_station_obs_table.py
"""

import glob
import gzip
import os

import numpy as np
import pandas as pd

BASE_DIR = "/private/groups/meyerlab/esteiger"
STATION_OBS_DIR = f"{BASE_DIR}/wind/station_obs"
OUT_OBS = f"{BASE_DIR}/wind/station_wind_obs.csv"

KT_TO_MS = 0.514444


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
    # dir 0 with speed 0 is calm, not north
    calm = (df["dir_deg"] == 0) & (df["speed_ms"] == 0)
    df.loc[calm, "dir_deg"] = np.nan
    return df[["station", "datetime", "dir_deg", "speed_ms"]]


def main():
    ptgc1 = parse_ptgc1()
    print(f"PTGC1: {len(ptgc1)} raw observations parsed")

    frames = [ptgc1]
    for st in ["VGN", "TQS"]:
        df = parse_asos(st)
        print(f"{st}: {len(df)} raw observations parsed")
        frames.append(df)

    obs = pd.concat(frames, ignore_index=True)
    obs.to_csv(OUT_OBS, index=False)
    print(f"Saved combined observation table (all months) -> {OUT_OBS} ({len(obs)} rows)")


if __name__ == "__main__":
    main()
