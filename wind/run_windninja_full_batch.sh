#!/bin/bash
# Run WindNinja terrain downscaling over the full 6-population study region
# for all 1464 bloom-season climatology timesteps, producing per-timestep
# speed+direction grids. Used to build a WindNinja-direction + HRRR-speed
# hybrid wind climatology (see PTGC1 validation: WindNinja improves direction
# but not speed at an exposed cape station).
set -euo pipefail

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate windninja

WORKDIR=/private/groups/meyerlab/esteiger/wind/windninja_full
DEM="${WORKDIR}/region_dem.tif"

run_one() {
  f="$1"
  WindNinja_cli \
    --num_threads 2 \
    --elevation_file "$DEM" \
    --initialization_method wxModelInitialization \
    --time_zone America/Los_Angeles \
    --forecast_filename "$f" \
    --mesh_resolution 460 --units_mesh_resolution m \
    --output_speed_units mps \
    --output_wind_height 10 --units_output_wind_height m \
    --vegetation grass \
    --write_ascii_output true \
    > /dev/null 2>>"${WORKDIR}/windninja_batch_errors.log"
}
export -f run_one
export DEM WORKDIR

find "${WORKDIR}/grib2" -name "*.grib2" | \
  xargs -P 16 -I{} bash -c 'run_one "$@"' _ {}

echo "Batch done. vel.asc outputs:"
find "${WORKDIR}/grib2" -name "*_vel.asc" | wc -l
