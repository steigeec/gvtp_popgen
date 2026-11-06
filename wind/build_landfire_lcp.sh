#!/bin/bash
# Build the LANDFIRE landscape (.lcp) files for WindNinja
# (run_windninja_landfire_batch.sh, run_windninja_tranquillon_sensitivity.sh).
#
# Layers on a 30 m grid, EPSG:32610:
#   elevation, slope, aspect  USGS 3DEP 1/3 arc-second DEM resampled to 30 m
#   fuel model                LANDFIRE LF2016 FBFM40
#   canopy                    LANDFIRE LF2016 CC (%), CH (m x 10), CBH (m x 10),
#                             CBD (kg/m3 x 100)
# LANDFIRE nodata (ocean) is set to fuel model 98 (open water) with no canopy.
#
# Output: BASE_DIR/env/landfire/
#   region_landfire.lcp       full study region, same extent as region_dem.tif
#   tranquillon_landfire.lcp  Tranquillon Mtn samples plus 2.5 km on each side
#
# LFPS needs an email with each request; pass it in LFPS_EMAIL.
# Env:    internet access (login node), GDAL (windninja conda env),
#         R with terra >= 1.8 (R_libs_windscape)
# Usage:  LFPS_EMAIL=you@example.edu bash build_landfire_lcp.sh
set -euo pipefail

export BASE_DIR=/private/groups/meyerlab/esteiger
: "${LFPS_EMAIL:?set LFPS_EMAIL to the email address LFPS requires}"

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate windninja

OUTDIR=${BASE_DIR}/env/landfire
mkdir -p "${OUTDIR}"
TMP=$(mktemp -d "${OUTDIR}/tmp_XXXX")
trap 'rm -rf "${TMP}"' EXIT

LAYERS="LF2016_FBFM40;LF2016_CC;LF2016_CH;LF2016_CBH;LF2016_CBD"
AOI="-120.77 34.37 -120.08 34.98"   # lon/lat, padded around the region domain
REGION_TE="706397.728179744 3809933.699465313 765637.728179744 3870203.699465313"
TRANQ_TE="723000 3825780 730770 3833550"
TILE=https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/13/TIFF/current/n35w121/USGS_13_n35w121.tif

# ---------- 1. LANDFIRE via LFPS ----------
job=$(curl -sS --max-time 120 -G https://lfps.usgs.gov/api/job/submit \
  --data-urlencode "Email=${LFPS_EMAIL}" --data-urlencode "Layer_List=${LAYERS}" \
  --data-urlencode "Area_of_Interest=${AOI}" --data-urlencode "Output_Projection=32610" \
  --data-urlencode "Resample_Resolution=30" | python3 -c "import json,sys; print(json.load(sys.stdin)['jobId'])")
echo "LFPS job ${job}"
while :; do
  status_json=$(curl -sS --max-time 60 "https://lfps.usgs.gov/api/job/status?JobId=${job}")
  status=$(echo "${status_json}" | python3 -c "import json,sys; print(json.load(sys.stdin)['status'])")
  [ "${status}" = "Succeeded" ] && break
  [ "${status}" = "Failed" ] && { echo "LFPS job failed"; exit 1; }
  sleep 15
done
url=$(echo "${status_json}" | python3 -c "import json,sys; print(json.load(sys.stdin)['outputFile'])")
curl -sS --max-time 600 -o "${TMP}/lf.zip" "${url}"
unzip -q -o "${TMP}/lf.zip" -d "${TMP}/lf"
LF_TIF=$(ls "${TMP}"/lf/*.tif)

# ---------- 2. Terrain on the LANDFIRE grid ----------
read -r xmin ymax xres <<<"$(gdalinfo -json "${LF_TIF}" | python3 -c "
import json,sys; g=json.load(sys.stdin)['geoTransform']; print(g[0], g[3], g[1])")"
read -r nx ny <<<"$(gdalinfo -json "${LF_TIF}" | python3 -c "
import json,sys; s=json.load(sys.stdin)['size']; print(s[0], s[1])")"
xmax=$(python3 -c "print(${xmin} + ${nx} * ${xres})"); ymin=$(python3 -c "print(${ymax} - ${ny} * ${xres})")
gdalwarp -q -overwrite -t_srs EPSG:32610 -te "${xmin}" "${ymin}" "${xmax}" "${ymax}" -tr 30 30 \
  -r bilinear "/vsicurl/${TILE}" "${TMP}/elev.tif"
gdaldem slope -q "${TMP}/elev.tif" "${TMP}/slope.tif"
gdaldem aspect -q -zero_for_flat "${TMP}/elev.tif" "${TMP}/aspect.tif"

# ---------- 3. Stack, fill ocean nodata ----------
Rscript - "${TMP}" "${LF_TIF}" <<'EOF'
.libPaths(file.path(Sys.getenv("BASE_DIR"), c("software/R_libs_windscape", "software/R_libs/")))
suppressMessages(library(terra))
a <- commandArgs(trailingOnly = TRUE)
lf <- rast(a[2])
fuel <- subst(lf[[1]], NA, 98)                  # ocean: fuel 98 (open water)
canopy <- subst(lf[[2:5]], NA, 0)               # and no canopy
topo <- c(rast(file.path(a[1], "elev.tif")), rast(file.path(a[1], "slope.tif")), rast(file.path(a[1], "aspect.tif")))
topo <- subst(round(topo), NA, 0)               # gdaldem leaves the outermost ring NA
stack <- c(topo, fuel, canopy)
names(stack) <- c("elevation", "slope", "aspect", "fuel", "cc", "ch", "cbh", "cbd")
writeRaster(stack, file.path(a[1], "stack.tif"), datatype = "INT2S", overwrite = TRUE, NAflag = -9999)
EOF

# ---------- 4. Write .lcp per domain ----------
LCP_OPTS="-co ELEVATION_UNIT=METERS -co SLOPE_UNIT=DEGREES -co ASPECT_UNIT=AZIMUTH_DEGREES \
  -co FUEL_MODEL_OPTION=NO_CUSTOM_AND_NO_FILE -co CANOPY_COV_UNIT=PERCENT \
  -co CANOPY_HT_UNIT=METERS_X_10 -co CBH_UNIT=METERS_X_10 -co CBD_UNIT=KG_PER_CUBIC_METER_X_100 \
  -co LATITUDE=34"
for d in "region ${REGION_TE}" "tranquillon ${TRANQ_TE}"; do
  set -- ${d}
  gdalwarp -q -overwrite -te "$2" "$3" "$4" "$5" -tr 30 30 -r near "${TMP}/stack.tif" "${TMP}/$1.tif"
  gdal_translate -q -of LCP ${LCP_OPTS} "${TMP}/$1.tif" "${OUTDIR}/$1_landfire.lcp"
  gdalinfo "${OUTDIR}/$1_landfire.lcp" | grep -E "Size is|Pixel Size"
done
echo "Done -> ${OUTDIR}"
