#!/bin/bash
# WindNinja sensitivity runs on the Tranquillon Mtn domain. Each variant
# changes one setting from the base run (460 m mesh, grass, no diurnal winds,
# neutral stability, 3DEP terrain):
#   base        base settings
#   mesh100     100 m mesh
#   landfire    tranquillon_landfire.lcp instead of grass
#   diurnal     --diurnal_winds true
#   stability   --non_neutral_stability true
#
# Input:  HRRR grib2 subsets (download_hrrr_full_subset.py), region_dem.tif
#         (build_windninja_dem.sh), tranquillon_landfire.lcp (build_landfire_lcp.sh)
# Output: BASE_DIR/wind/windninja_landfire/sensitivity/<variant>/,
#         read by windscape_tranquillon_sensitivity.R
# Env:    windninja conda env
# Usage:  sbatch --array=0-4 run_windninja_tranquillon_sensitivity.sh
# Timesteps already run are skipped.

#SBATCH --job-name=wn_sens
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=4:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/wind/windninja_landfire/sensitivity/wn_sens_%A_%a.log

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

VARIANTS=(base mesh100 landfire diurnal stability)
VARIANT="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate windninja

WIND_DIR=${BASE_DIR}/wind
GRIB_DIR="${WIND_DIR}/windninja_full/grib2"
SENS_DIR="${WIND_DIR}/windninja_landfire/sensitivity"
OUT_DIR="${SENS_DIR}/${VARIANT}"
mkdir -p "${OUT_DIR}"

# Tranquillon crop of the 3DEP DEM, used by all variants except landfire
DEM="${SENS_DIR}/tranquillon_dem.tif"
if [ ! -s "${DEM}" ]; then
  gdal_translate -q -projwin 723000 3833550 730770 3825780 \
    ${BASE_DIR}/env/DEM/region_dem.tif "${DEM}.${VARIANT}.tmp.tif" && mv "${DEM}.${VARIANT}.tmp.tif" "${DEM}"
fi

ELEV="${DEM}"; VEG="--vegetation grass"; MESH=460; DIURNAL=false; STABILITY=false
case "${VARIANT}" in
  base) ;;
  mesh100)   MESH=100 ;;
  landfire)  ELEV=${BASE_DIR}/env/landfire/tranquillon_landfire.lcp; VEG="" ;;
  diurnal)   DIURNAL=true ;;
  stability) STABILITY=true ;;
esac

run_one() {
  f="$1"
  ymd=$(basename "$(dirname "$f")")
  hh=$(basename "$f" | sed -E 's/.*hrrr\.t([0-9]{2})z.*/\1/')
  local_t=$(TZ=America/Los_Angeles date -d "${ymd:0:4}-${ymd:4:2}-${ymd:6:2} ${hh}:00 UTC" +%m-%d-%Y_%H%M)
  ls "${OUT_DIR}"/*_"${local_t}"_*vel.asc >/dev/null 2>&1 && return 0
  WindNinja_cli \
    --num_threads 2 \
    --elevation_file "${ELEV}" ${VEG} \
    --initialization_method wxModelInitialization \
    --time_zone America/Los_Angeles \
    --forecast_filename "$f" \
    --mesh_resolution "${MESH}" --units_mesh_resolution m \
    --diurnal_winds "${DIURNAL}" \
    --non_neutral_stability "${STABILITY}" \
    --output_speed_units mps \
    --output_wind_height 10 --units_output_wind_height m \
    --write_ascii_output true \
    --output_path "${OUT_DIR}" \
    > /dev/null 2>>"${OUT_DIR}/windninja_errors.log"
}
export -f run_one
export ELEV VEG MESH DIURNAL STABILITY OUT_DIR

find "${GRIB_DIR}" -name "*.grib2" | xargs -P 16 -I{} bash -c 'run_one "$@"' _ {}

echo "Batch done. ${VARIANT} vel.asc outputs: $(find "${OUT_DIR}" -name '*_vel.asc' | wc -l)"
