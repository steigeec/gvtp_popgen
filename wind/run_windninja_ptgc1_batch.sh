#!/bin/bash
# Run WindNinja terrain downscaling for all 486 HRRR timesteps matched
# against PTGC1 (Point Arguello) observations, to test whether downscaling
# closes the ~28% raw-HRRR speed underestimate found there. Each timestep's
# 4-band grib2 subset already downloaded by download_hrrr_ptgc1_subset.py.
set -euo pipefail

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate windninja

WORKDIR=/private/groups/meyerlab/esteiger/wind/windninja_ptgc1
DEM="${WORKDIR}/ptgc1_dem.tif"

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
