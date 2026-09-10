#!/bin/bash
# Download raw wind direction/speed observations for stations bracketing
# the Point Conception cape, for building a directional wind climatology
# (wind rose) per population site. Run on a node with internet access
# (not a compute node without egress) -- this is a one-time interactive
# download, not a SLURM job.
#
# Station network (chosen as nearest-station per site; see
# gvtp_popgen/wind/wind_connectivity.R header for the full distance table):
#   PTGC1 - NDBC C-MAN station AT Point Arguello (34.577,-120.648).
#           Nearest station for VAFB_Arguello (1.5 km), Dangermond (13.9 km),
#           Point_Conception (18.7 km). Land station, anemometer 7 m AGL.
#           Discontinued 2023-02-14; historical record used through then.
#   VGN   - Vandenberg SFB North ASOS (34.8495,-120.5996).
#           Nearest station for VAFB_Lion_Head (1.6 km).
#   TQS   - Vandenberg SFB South ASOS (34.636,-120.6161).
#           Nearest/tied-nearest for Strauss (10.6 km).
#   SBA   - Santa Barbara Municipal ASOS (34.4261,-119.8436).
#           Best available proxy for CDFW (32.7 km) -- CDFW sits on the lee
#           side of the cape in the Santa Barbara Channel, and SBA is the
#           closest coastal station sampling that same sheltered regime.
#           No closer station exists; this is a known limitation.
#
# Source: NDBC historical standard meteorological data (public, no auth)
# and Iowa Environmental Mesonet ASOS archive (public, no auth).

set -euo pipefail

OUTDIR=/private/groups/meyerlab/esteiger/wind/station_obs
mkdir -p "${OUTDIR}"

# ---------- NDBC: PTGC1 historical standard meteorological data ----------
# Yearly gzipped files; not every year exists (station commissioned partway
# through the archive, decommissioned 2023-02-14) so 404s are expected and
# skipped rather than treated as errors.
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

# ---------- IEM ASOS: VGN, TQS, SBA hourly direction + speed ----------
echo "=== IEM ASOS: VGN, TQS, SBA (2000-2024) ==="
for station in VGN TQS SBA; do
  out="${OUTDIR}/${station}_asos_2000_2024.csv"
  url="https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=${station}&data=drct,sknt&year1=2000&month1=1&day1=1&year2=2024&month2=12&day2=31&tz=UTC&format=onlycomma&latlon=no&elev=no&missing=empty&trace=empty&direct=no&report_type=3"
  curl -s "${url}" -o "${out}"
  n=$(wc -l < "${out}")
  echo "  ${station}: ${n} lines -> ${out}"
done

echo "Done. Raw station data in ${OUTDIR}"
