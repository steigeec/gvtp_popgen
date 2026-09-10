#!/private/home/esteiger/miniconda3/envs/herbie/bin/python3
"""
Batch-download the 4-band HRRR subset (U/V 10m, TMP 2m, TCDC entire atm --
the set WindNinja's HRRR reader requires) for ALL 1464 bloom-season
climatology timesteps (2022-2024, Jun-Sep, 6-hourly -- same set as
download_hrrr_wind.py), for the full-region WindNinja production run used to
get per-population WindNinja direction (see PTGC1 validation: WindNinja
direction is more accurate than raw HRRR direction there, so this is used to
build a WindNinja-direction + HRRR-speed hybrid wind climatology for
wind_connectivity.R / wind_aem_rda.R).
"""
import pandas as pd
from herbie import FastHerbie

YEARS = [2022, 2023, 2024]
BLOOM_MONTHS = [6, 7, 8, 9]
HOURS = [0, 6, 12, 18]
SAVE_DIR = "/private/groups/meyerlab/esteiger/wind/windninja_full/grib2"
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
