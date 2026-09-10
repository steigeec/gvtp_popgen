#!/private/home/esteiger/miniconda3/envs/herbie/bin/python3
"""
Build a per-population wind climatology that uses WindNinja's downscaled
DIRECTION (more accurate at exposed points per the PTGC1 validation: median
direction error 16 deg -> 10 deg) but raw HRRR SPEED (WindNinja doesn't
recover HRRR's speed underestimate there -- 28% -> 32% low, i.e. no better,
so there's no reason to trade HRRR's better-validated speed for WindNinja's).

Unlike process_hrrr_climatology.py (per HRRR-grid-cell, 478 cells), this is
built directly per POPULATION (6 sites): each population samples WindNinja's
own output grid at its own coordinates (not a shared grid cell), and pairs
that with HRRR speed at that population's own nearest HRRR grid cell, for
every one of the 1464 bloom-season timesteps that has a WindNinja run.

Output mirrors hrrr_windrose_summary.txt / hrrr_windrose_sectors.txt's
schema (same column names) but with cell_id = population code, so
wind_connectivity.R / wind_aem_rda.R need only a trivial pop->cell_id change
(identity instead of nearest_cell()) to consume it -- see
hybrid_wind_connectivity.R.
"""
import glob
import re
import subprocess

import numpy as np
import pandas as pd

WIND_DIR = "/private/groups/meyerlab/esteiger/wind"
WN_DIR = f"{WIND_DIR}/windninja_full/grib2"
HRRR_RAW_CSV = f"{WIND_DIR}/hrrr/hrrr_wind_raw_points.csv"
OUT_SUMMARY = f"{WIND_DIR}/hrrr/hybrid_windrose_summary.txt"
OUT_SECTORS = f"{WIND_DIR}/hrrr/hybrid_windrose_sectors.txt"

N_SECTORS = 16
SECTOR_WIDTH = 360 / N_SECTORS

NAME_MAP = {
    "Point_Conception": "conception", "CDFW": "cdfw",
    "VAFB_Arguello": "arguello", "Strauss": "strauss",
    "VAFB_Lion_Head": "lionhead", "Dangermond": "Dangermond",
}


def to_from_bearing_uv(u, v):
    dir_to = np.degrees(np.arctan2(u, v)) % 360
    return (dir_to + 180) % 360


def load_population_centroids():
    field = pd.read_csv("/private/groups/meyerlab/esteiger/field_data.txt", sep="\t")
    cent = field.groupby("Clusters")[["Lat", "Long"]].mean().reset_index()
    cent["pop"] = cent["Clusters"].map(NAME_MAP)
    return cent


def load_hrrr_raw():
    df = pd.read_csv(HRRR_RAW_CSV, parse_dates=["datetime"])
    df["lon180"] = df["lon"] - 360
    return df


def nearest_hrrr_cell(hrrr, lat, lon):
    first_t = hrrr["datetime"].iloc[0]
    snap = hrrr[hrrr["datetime"] == first_t]
    d = np.hypot(snap["lat"] - lat, snap["lon180"] - lon)
    cell_lat, cell_lon = snap.loc[d.idxmin(), ["lat", "lon180"]]
    sub = hrrr[(hrrr["lat"] == cell_lat) & (hrrr["lon180"] == cell_lon)].copy()
    sub["speed_hrrr"] = np.hypot(sub["u10"], sub["v10"])
    return sub[["datetime", "speed_hrrr"]]


def sample_asc(path, lat, lon):
    out = subprocess.run(
        ["gdallocationinfo", "-wgs84", "-valonly", path, str(lon), str(lat)],
        capture_output=True, text=True,
    )
    v = out.stdout.strip()
    return float(v) if v else np.nan


def parse_dt_from_name(path):
    m = re.search(r"_(\d{2})-(\d{2})-(\d{4})_(\d{2})(\d{2})_460m_ang\.asc$", path)
    mo, d, y, hh, mm = m.groups()
    local = pd.Timestamp(int(y), int(mo), int(d), int(hh), int(mm), tz="America/Los_Angeles")
    return local.tz_convert("UTC").tz_localize(None)


def main():
    ang_files = sorted(glob.glob(f"{WN_DIR}/**/*_ang.asc", recursive=True))
    print(f"Found {len(ang_files)} WindNinja direction grids")
    wn_times = pd.DataFrame({"datetime": [parse_dt_from_name(f) for f in ang_files],
                              "ang_path": ang_files})

    centroids = load_population_centroids()
    hrrr = load_hrrr_raw()

    summary_rows, sector_frames = [], []
    for _, row in centroids.iterrows():
        pop, lat, lon = row["pop"], row["Lat"], row["Long"]
        hrrr_speed = nearest_hrrr_cell(hrrr, lat, lon)
        merged = wn_times.merge(hrrr_speed, on="datetime", how="inner")
        merged["dir_wn"] = [sample_asc(p, lat, lon) for p in merged["ang_path"]]
        merged = merged.dropna(subset=["dir_wn"])

        speed = merged["speed_hrrr"].values
        dir_from = merged["dir_wn"].values
        u = -speed * np.sin(np.radians(dir_from))
        v = -speed * np.cos(np.radians(dir_from))
        mean_speed = speed.mean()
        mean_u, mean_v = u.mean(), v.mean()
        resultant_speed = np.hypot(mean_u, mean_v)
        R = resultant_speed / mean_speed if mean_speed > 0 else np.nan
        resultant_dir = to_from_bearing_uv(np.array([mean_u]), np.array([mean_v]))[0]

        summary_rows.append((pop, lat, lon, len(merged), round(mean_speed, 3),
                              round(resultant_dir, 1), round(R, 3)))

        sector = np.floor((dir_from + SECTOR_WIDTH / 2) / SECTOR_WIDTH).astype(int) % N_SECTORS
        for s in range(N_SECTORS):
            in_sector = speed[sector == s]
            freq = len(in_sector) / len(merged)
            mean_spd_sector = in_sector.mean() if len(in_sector) else 0.0
            sector_frames.append((pop, s * SECTOR_WIDTH, freq, mean_spd_sector, len(in_sector)))

        print(f"  {pop}: {len(merged)} matched timesteps, "
              f"mean speed (HRRR) {mean_speed:.2f} m/s, "
              f"resultant dir (WindNinja) {resultant_dir:.0f} deg, R={R:.3f}")

    summary = pd.DataFrame(summary_rows, columns=[
        "cell_id", "lat", "lon", "n_obs", "mean_speed_ms",
        "resultant_dir_deg_from", "resultant_R"])
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
    print(f"\nSaved hybrid per-population summary -> {OUT_SUMMARY}")

    sectors = pd.DataFrame(sector_frames, columns=[
        "cell_id", "sector_center_deg", "frequency", "mean_speed_ms", "n_obs"])
    sectors.to_csv(OUT_SECTORS, sep="\t", index=False)
    print(f"Saved hybrid per-population sector tables -> {OUT_SECTORS}")


if __name__ == "__main__":
    main()
