#!/private/home/esteiger/miniconda3/envs/herbie/bin/python3
"""
Batch-download the 4-band HRRR subset (U/V 10m, TMP 2m, TCDC entire atm --
the set WindNinja's HRRR reader requires, see
gvtp_popgen/wind/wind[...]/herbie gotchas) for the exact 486 timesteps
matched against PTGC1 obs in hrrr/validate_PTGC1_matched.csv (2022
bloom-season, 6-hourly). Used to test whether WindNinja downscaling closes
the ~28% raw-HRRR speed underestimate found at PTGC1 (Point Arguello).
"""
import pandas as pd
from herbie import FastHerbie

MATCHED_CSV = "/private/groups/meyerlab/esteiger/wind/hrrr/validate_PTGC1_matched.csv"
SAVE_DIR = "/private/groups/meyerlab/esteiger/wind/windninja_ptgc1/grib2"
BANDS = ":(UGRD|VGRD):10 m above ground|:TMP:2 m above ground|:TCDC:entire atmosphere"

times = pd.read_csv(MATCHED_CSV, parse_dates=["datetime"])["datetime"].tolist()
print(f"Downloading {len(times)} timesteps...")

FH = FastHerbie(times, model="hrrr", product="sfc", fxx=[0], max_threads=30)
paths = FH.download(BANDS, max_threads=30, save_dir=SAVE_DIR)
print(f"Downloaded {len(paths)} files -> {SAVE_DIR}")
