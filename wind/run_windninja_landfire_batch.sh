#!/bin/bash
# WindNinja downscaling of HRRR for all 4896 bloom-season timesteps on a
# LANDFIRE landscape (build_landfire_lcp.sh), with diurnal winds and
# non-neutral stability.
#
# Input:  HRRR grib2 subsets (download_hrrr_full_subset.py),
#         BASE_DIR/env/landfire/<domain>_landfire.lcp
# Output: BASE_DIR/wind/windninja_landfire/<domain>/*_vel.asc, *_ang.asc
# Env:    windninja conda env
# Usage:  sbatch run_windninja_landfire_batch.sh region       (460 m mesh)
#         sbatch run_windninja_landfire_batch.sh tranquillon  (100 m mesh)
# Timesteps already run are skipped.

#SBATCH --job-name=windninja_lf
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/wind/windninja_landfire_%x_%j.log

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

DOMAIN="${1:?usage: run_windninja_landfire_batch.sh region|tranquillon}"
case "${DOMAIN}" in
  region)      MESH=460 ;;
  tranquillon) MESH=100 ;;
  *) echo "unknown domain ${DOMAIN}"; exit 1 ;;
esac

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate windninja

WIND_DIR=${BASE_DIR}/wind
GRIB_DIR="${WIND_DIR}/windninja_full/grib2"
LCP="${BASE_DIR}/env/landfire/${DOMAIN}_landfire.lcp"
OUT_DIR="${WIND_DIR}/windninja_landfire/${DOMAIN}"
mkdir -p "${OUT_DIR}"

# WindNinja output name, used by the skip check: <lcp basename>_<local MM-DD-YYYY_HHMM>_<mesh>m_non_neutral_stability_*.asc
# e.g. grib2/hrrr/20220601/...hrrr.t00z... gives region_landfire_05-31-2022_1700_460m_non_neutral_stability_vel.asc
run_one() {
  f="$1"
  ymd=$(basename "$(dirname "$f")")
  hh=$(basename "$f" | sed -E 's/.*hrrr\.t([0-9]{2})z.*/\1/')
  local_t=$(TZ=America/Los_Angeles date -d "${ymd:0:4}-${ymd:4:2}-${ymd:6:2} ${hh}:00 UTC" +%m-%d-%Y_%H%M)
  [ -s "${OUT_DIR}/${DOMAIN}_landfire_${local_t}_${MESH}m_non_neutral_stability_vel.asc" ] && return 0
  WindNinja_cli \
    --num_threads 2 \
    --elevation_file "${LCP}" \
    --initialization_method wxModelInitialization \
    --time_zone America/Los_Angeles \
    --forecast_filename "$f" \
    --mesh_resolution "${MESH}" --units_mesh_resolution m \
    --diurnal_winds true \
    --non_neutral_stability true \
    --output_speed_units mps \
    --output_wind_height 10 --units_output_wind_height m \
    --write_ascii_output true \
    --output_path "${OUT_DIR}" \
    > /dev/null 2>>"${OUT_DIR}/windninja_errors.log"
}
export -f run_one
export DOMAIN MESH LCP OUT_DIR

find "${GRIB_DIR}" -name "*.grib2" | xargs -P 16 -I{} bash -c 'run_one "$@"' _ {}

echo "Batch done. ${DOMAIN} vel.asc outputs:"
find "${OUT_DIR}" -name "*_vel.asc" | wc -l
