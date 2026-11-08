#!/bin/bash
# Fetch Deinandra occurrence records (GBIF API, curl + jq), a Natural Earth
# land polygon, and the manual Jepson Herbarium (CCH2) export, then run
# phylogeny_range.R to draw the phylogeny + range figure.
# GBIF fetch needs internet access; run once on the login node to cache
# occurrences.csv if compute nodes are offline.
#
# Usage: bash phylogeny-range.sh

#SBATCH --job-name=phylogeny_range
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=0:30:00
#SBATCH --output=phylogeny_range_%j.log

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

SCRIPT_DIR=${BASE_DIR}/gvtp_popgen/phylogeny_range
OUTDIR=${BASE_DIR}/phylogeny_range
mkdir -p "$OUTDIR"

TREE_FILE="$OUTDIR/deinandra_phylogeny.nwk"
OCC_CSV="$OUTDIR/deinandra_occurrences.csv"
VILLOSA_CSV=${BASE_DIR}/pheno/updatedvarsandpheno.csv
OUT_PNG="$OUTDIR/deinandra_phylogeny_range.png"
BASEMAP_DIR="$OUTDIR/basemap"
LAND_SHP="$BASEMAP_DIR/ne_10m_land.shp"

# Manual CCH2 export for UCJEPS (not mirrored in GBIF; CCH2 has no query API),
# downloaded from https://www.cch2.org/portal/
JEPSON_ZIP="$OUTDIR/SymbOutput_2026-09-05_194011_DwC-A.zip"
JEPSON_DIR="$OUTDIR/jepson_dwca"
JEPSON_CSV="$JEPSON_DIR/occurrences.csv"

# GBIF genus key for Deinandra
GENUS_KEY=3106606

# Map query window
LAT_RANGE="33.9,35.3"
LON_RANGE="-121.3,-119.5"

# --- Fetch Deinandra occurrence records from GBIF (cached; skip if present) --
if [ -f "$OCC_CSV" ]; then
  echo "$(date) -- Occurrence cache exists, skipping GBIF fetch: $OCC_CSV"
else
  echo "$(date) -- Fetching Deinandra occurrences from GBIF ..."
  echo "species,infraspecificEpithet,decimalLatitude,decimalLongitude,coordinateUncertaintyInMeters,institutionCode" > "$OCC_CSV"

  OFFSET=0
  LIMIT=300
  while true; do
    RESP=$(curl -s "https://api.gbif.org/v1/occurrence/search?taxonKey=${GENUS_KEY}&hasCoordinate=true&decimalLatitude=${LAT_RANGE}&decimalLongitude=${LON_RANGE}&limit=${LIMIT}&offset=${OFFSET}")

    N=$(echo "$RESP" | jq '.results | length')
    echo "$RESP" | jq -r '.results[] |
      [(.species // ""), (.infraspecificEpithet // ""), (.decimalLatitude // ""), (.decimalLongitude // ""),
       (.coordinateUncertaintyInMeters // ""), (.institutionCode // "")]
      | @csv' >> "$OCC_CSV"

    echo "$(date) --   offset ${OFFSET}: ${N} records"
    if [ "$N" -lt "$LIMIT" ]; then
      break
    fi
    OFFSET=$((OFFSET + LIMIT))
  done
  echo "$(date) -- Done fetching GBIF occurrences: $(wc -l < "$OCC_CSV") lines"
fi

# --- Fetch a high-resolution coastline (cached) ------------------------------
# the maps package coastline is too coarse at Point Conception
if [ -f "$LAND_SHP" ]; then
  echo "$(date) -- Basemap cache exists, skipping download: $LAND_SHP"
else
  echo "$(date) -- Fetching Natural Earth 10m land polygon ..."
  mkdir -p "$BASEMAP_DIR"
  curl -sL -o "$BASEMAP_DIR/ne_10m_land.zip" "https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_land.zip"
  unzip -o -q "$BASEMAP_DIR/ne_10m_land.zip" -d "$BASEMAP_DIR"
fi

# --- Unpack the manual Jepson Herbarium (CCH2) export (cached) --------------
if [ -f "$JEPSON_CSV" ]; then
  echo "$(date) -- Jepson DwC-A cache exists, skipping unzip: $JEPSON_CSV"
else
  echo "$(date) -- Unpacking Jepson Herbarium DwC-A export ..."
  mkdir -p "$JEPSON_DIR"
  unzip -o -q "$JEPSON_ZIP" -d "$JEPSON_DIR"
fi

# --- Build the figure --------------------------------------------------------
echo "$(date) -- Rendering phylogeny + range figure ..."
Rscript "$SCRIPT_DIR/phylogeny_range.R" "$TREE_FILE" "$OCC_CSV" "$VILLOSA_CSV" "$JEPSON_CSV" "$LAND_SHP" "$OUT_PNG"
echo "$(date) -- Done: $OUT_PNG"
