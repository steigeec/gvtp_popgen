#!/bin/bash
# Download the OpenStreetMap coastline (natural=coastline) for the study area
# from the Overpass API, for the wind maps (windscape_composite_figure.R).
#
# Output: BASE_DIR/env/coastline/osm_coastline_study_area.gpkg
# Env:    internet access (login node), GDAL (windninja conda env)
# Usage:  bash build_coastline.sh
set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate windninja

OUTDIR=${BASE_DIR}/env/coastline
mkdir -p "${OUTDIR}"
BBOX="34.2,-121.0,35.1,-119.8"   # south, west, north, east

TMP=$(mktemp -d "${OUTDIR}/tmp_osm_XXXX")
trap 'rm -rf "${TMP}"' EXIT
curl -sS --fail --max-time 300 -A "gvtp_popgen coastline (curl)" -o "${TMP}/coastline.osm" \
  --data-urlencode "data=[out:xml][timeout:120];way[\"natural\"=\"coastline\"](${BBOX});(._;>;);out body;" \
  https://overpass-api.de/api/interpreter

ogr2ogr -overwrite -f GPKG "${OUTDIR}/osm_coastline_study_area.gpkg" \
  "${TMP}/coastline.osm" lines -nln coastline

ogrinfo -so "${OUTDIR}/osm_coastline_study_area.gpkg" coastline | grep -E "Feature Count|Extent"
echo "Done -> ${OUTDIR}/osm_coastline_study_area.gpkg"
