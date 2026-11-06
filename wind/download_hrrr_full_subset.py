#!/private/home/esteiger/miniconda3/envs/herbie/bin/python3
"""
Download the 4-band HRRR subset WindNinja needs (UGRD/VGRD 10 m, TMP 2 m,
TCDC entire atmosphere) for the same 4896 timesteps as download_hrrr_wind.py.

Output: GRIB2 files in BASE_DIR/wind/windninja_full/grib2,
        read by run_windninja_landfire_batch.sh and
        run_windninja_tranquillon_sensitivity.sh
Env:    'herbie' conda env
Usage:  /private/home/esteiger/miniconda3/envs/herbie/bin/python3 download_hrrr_full_subset.py
Files already on disk are skipped by Herbie.
"""
import pandas as pd
from herbie import FastHerbie

YEARS = list(range(2015, 2023))   # pre-turbine (Strauss Wind Farm began operating 2023)
BLOOM_MONTHS = [6, 7, 8, 9, 10]   # Gaviota tarplant flowering season, June-October
HOURS = [0, 6, 12, 18]
BASE_DIR = "/private/groups/meyerlab/esteiger"
SAVE_DIR = f"{BASE_DIR}/wind/windninja_full/grib2"
BANDS = ":(UGRD|VGRD):10 m above ground|:TMP:2 m above ground|:TCDC:entire atmosphere"


def build_datetimes():
    out = []
    for year in YEARS:
        for month in BLOOM_MONTHS:
            days = pd.Period(f"{year}-{month:02d}").days_in_month
            out += [pd.Timestamp(year=year, month=month, day=day, hour=hour)
                    for day in range(1, days + 1) for hour in HOURS]
    return out


times = build_datetimes()
print(f"Downloading {len(times)} timesteps...")

FH = FastHerbie(times, model="hrrr", product="sfc", fxx=[0], max_threads=40)
paths = FH.download(BANDS, max_threads=40, save_dir=SAVE_DIR)
print(f"Downloaded {len(paths)} files -> {SAVE_DIR}")
