#!/bin/bash
# Download raw wind observations for the HRRR and WindNinja validation
# (validate_hrrr_vs_stations.py, validate_windninja_stations.py).
#
# Stations:
#   PTGC1  NDBC C-MAN, Point Arguello (34.577, -120.648), anemometer 7 m AGL;
#          discontinued 2023-02-14
#   VGN    Vandenberg SFB North ASOS (34.8495, -120.5996)
#   TQS    Vandenberg SFB South ASOS (34.636, -120.6161)
# Sources: NDBC historical stdmet files and the Iowa Environmental Mesonet
# ASOS archive (both public).
#
# Output: BASE_DIR/wind/station_obs/
# Env:    internet access; run interactively, not as a SLURM job
# Usage:  bash download_station_wind.sh

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

OUTDIR=${BASE_DIR}/wind/station_obs
mkdir -p "${OUTDIR}"

# ---------- NDBC: PTGC1 ----------
# One file per year; years with no file return 404 and are skipped.
echo "=== NDBC PTGC1 (Point Arguello C-MAN) ==="
mkdir -p "${OUTDIR}/ptgc1"
for year in $(seq 1995 2023); do
  url="https://www.ndbc.noaa.gov/data/historical/stdmet/ptgc1h${year}.txt.gz"
  out="${OUTDIR}/ptgc1/ptgc1h${year}.txt.gz"
  code=$(curl -s -o "${out}" -w "%{http_code}" "${url}")
  if [ "${code}" = "200" ]; then
    echo "  ${year}: OK"
  else
    echo "  ${year}: not available (HTTP ${code}), skipping"
    rm -f "${out}"
  fi
done

# ---------- IEM ASOS: VGN, TQS ----------
echo "=== IEM ASOS: VGN, TQS (2000-2024) ==="
for station in VGN TQS; do
  out="${OUTDIR}/${station}_asos_2000_2024.csv"
  url="https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=${station}&data=drct,sknt&year1=2000&month1=1&day1=1&year2=2024&month2=12&day2=31&tz=UTC&format=onlycomma&latlon=no&elev=no&missing=empty&trace=empty&direct=no&report_type=3"
  curl -s "${url}" -o "${out}"
  n=$(wc -l < "${out}")
  echo "  ${station}: ${n} lines -> ${out}"
done

echo "Done. Raw station data in ${OUTDIR}"
