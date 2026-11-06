# gvtp_popgen
Structure and selection in Deinandra increscens villosa

## Wind connectivity and gene flow

Tests whether wind connectivity explains genetic differentiation and diversity
among populations (HRRR wind) and among individuals on Tranquillon Mountain
(WindNinja-downscaled wind), using windscape (Kling & Ackerly 2020, 2021), and
whether gene flow estimated with FRAME is biased downwind. Wind is a pre-turbine
climatology: June-October 2015-2022, 6-hourly (the Strauss Wind Farm began
operating in 2023).

Each script sets `BASE_DIR` (the project data directory) at the top; inputs and
outputs are read from and written to paths under it, not to this repository.

### Run order

Wind data (`wind/`)
1. `download_station_wind.sh`, `build_station_obs_table.py`: station observations (PTGC1, VGN, TQS)
2. `download_hrrr_wind.py`, then `rotate_hrrr_wind.py`: HRRR 10 m wind, rotated to true north
3. `validate_hrrr_vs_stations.py`: HRRR vs station observations
4. `build_windninja_dem.sh`, `build_coastline.sh`, `build_landfire_lcp.sh`: terrain, coastline and LANDFIRE 2016 landscape files (`build_landfire_lcp.sh` needs `LFPS_EMAIL` set)
5. `download_hrrr_full_subset.py`: HRRR grib2 files used to initialize WindNinja
6. `sbatch run_windninja_landfire_batch.sh region` and `sbatch run_windninja_landfire_batch.sh tranquillon`; `sbatch --array=0-4 run_windninja_tranquillon_sensitivity.sh`
7. `validate_windninja_stations.py`: WindNinja vs station observations

Connectivity and tests (`wind/`)

8. `windscape_build_series.R`: wind time series rasters for windscape
9. `windscape_connectivity.R`: wind roses and pairwise travel times
10. `windscape_tests.R`: isolation and diversity tests (populations and Tranquillon Mtn)
11. `windscape_population_random_walk.R`: random-walk connectivity at the population level
12. `Rscript windscape_tranquillon_sensitivity.R base mesh100 landfire diurnal stability full`: WindNinja configuration sensitivity
13. `windscape_composite_figure.R`: connectivity figure

Population Fst checks (`wind/`, `structure/`)
- `fst_downsample_launch.sh random|covmatch` (calls `fst_downsample_saf.sh`, `fst_downsample_pair.sh`): Fst after downsampling to test coverage effects
- `structure/mantel_IBD.R`, `structure/plot_Fst_heatmap.R`: isolation by distance and Fst heatmap

Directional gene flow (`migration/`)
- `create_grid.sh` (`create_grid.R`): spatial grid shared by FEEMS and FRAME
- `FEEMS.sh`, `FEEMS_fit.sh`, `FEEMS_year2.sh`, `FEEMS_year3.sh`: FEEMS by year
- `FRAME_fit.sh`, `FRAME_lambda_sensitivity.sh`: FRAME pooled across years
- `frame_wind_asymmetry.R`: gene flow asymmetry vs wind flow asymmetry, with a rotated-wind null (run after step 8)

### Environments
- `herbie` conda env (Herbie, pandas, scipy): HRRR download, rotation and validation
- `windninja` conda env (WindNinja_cli, GDAL): DEM, coastline, LANDFIRE and WindNinja runs; `validate_windninja_stations.py` needs it on `PATH`
- `feems_e` and `frame_e` conda envs: FEEMS and FRAME (run with `PYTHONNOUSERSITE=1`)
- R: windscape 2.0.0 and terra >= 1.8 (installed in `software/R_libs_windscape`), plus geosphere, vegan, sf, jpeg, ggplot2, data.table and dggridR
