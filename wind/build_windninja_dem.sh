#!/bin/bash
# Build the 10 m DEM for the WindNinja terrain-only sensitivity run
# (run_windninja_tranquillon_sensitivity.sh).
#
# Input:  USGS 3DEP 1/3 arc-second tile n35w121, read over HTTP with /vsicurl/
#         https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/current/
# Output: BASE_DIR/env/DEM/region_dem.tif
#         (EPSG:32610, 10 m, Float32, m a.s.l.)
# The extent matches the existing WindNinja grid so outputs stay aligned.
# WindNinja names its outputs after the DEM (region_dem_<date>_..._vel.asc) and
# downstream scripts parse that name, so keep the basename.
#
# Conda env: windninja
#   conda create -n windninja -c conda-forge windninja=3.13.0.1
#   (also provides GDAL: gdalwarp, gdal_translate, gdallocationinfo)
# Needs internet access (login node).
# Usage:  bash build_windninja_dem.sh
set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate windninja

OUTDIR=${BASE_DIR}/env/DEM
mkdir -p "${OUTDIR}"

TILE_URL=https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/current
TILES=(n35w121)   # covers 34-35 N, 120-121 W

SRC_VRT="${OUTDIR}/usgs_3dep_13_tiles.vrt"
gdalbuildvrt -overwrite "${SRC_VRT}" \
  $(for t in "${TILES[@]}"; do echo "/vsicurl/${TILE_URL}/${t}/USGS_13_${t}.tif"; done)

build_dem() {
  local out="$1"; shift
  echo "=== ${out} ==="
  gdalwarp -overwrite \
    -t_srs EPSG:32610 -te "$@" -tr 10 10 -r bilinear \
    -ot Float32 -dstnodata nan \
    -co COMPRESS=DEFLATE -co PREDICTOR=3 -co TILED=YES \
    "${SRC_VRT}" "${OUTDIR}/${out}"
  gdalinfo -stats "${OUTDIR}/${out}" | grep -E "Size is|Pixel Size|STATISTICS_(MINIMUM|MAXIMUM|VALID_PERCENT)"
  rm -f "${OUTDIR}/${out}.aux.xml"
}

build_dem region_dem.tif 706397.728179744 3809933.699465313 765637.728179744 3870203.699465313

rm -f "${SRC_VRT}"
echo "Done. DEM in ${OUTDIR}"
