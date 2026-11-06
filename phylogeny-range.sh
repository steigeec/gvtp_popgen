#!/bin/bash
# =============================================================================
# phylogeny-range.sh
#
# Builds a two-panel figure for Deinandra increscens subsp. villosa (Gaviota
# tarplant):
#   1. its phylogenetic placement among close Deinandra relatives, reproduced
#      (topology only) from a published figure -- see $OUTDIR/deinandra_phylogeny.nwk
#      and phylogeny_range.R headers for the full citation trail (Simpson 2010
#      Plant Systematics Fig. 19.9, drawing on Baldwin's molecular work on
#      Deinandra);
#   2. a range map, zoomed out from (but centered on) this project's study
#      bounding box (bbox_extent <- ext(-120.7, -120.1, 34.4, 34.9); see e.g.
#      pheno/fire.R), showing GBIF occurrence records for D. increscens subsp.
#      villosa and any other Deinandra species/subspecies recorded in the
#      surrounding area.
#
# Occurrence data comes live from the GBIF API (taxonKey 3106606 = genus
# Deinandra) -- no rgbif/geodata R packages required, just curl + jq. The
# fetch step needs outbound internet access; if a compute node lacks it, run
# this script once on the login node so occurrences.csv gets cached, then
# reruns (including array/batch resubmits) will skip straight to plotting.
#
# D. increscens subsp. villosa is federally endangered, and iNaturalist
# automatically obscures coordinates for threatened taxa: essentially every
# iNaturalist record of it carries coordinateUncertaintyInMeters ~28800 (a
# ~29 km randomization radius, per each record's own informationWithheld
# field, e.g. "Coordinate uncertainty increased to 28818m to protect
# threatened taxon"). Herbarium specimen records (SBBG, UCSB, etc.) are not
# obscured this way, but between the two, GBIF's usable villosa sample is
# thin -- so villosa's points come from this project's own field/herbarium
# dataset instead (VILLOSA_CSV below), and GBIF only supplies the other
# Deinandra taxa. phylogeny_range.R still drops any GBIF record with >1 km
# coordinate uncertainty before plotting, so fetch coordinateUncertaintyInMeters
# and institutionCode for that filter.
#
# GBIF has near-complete coverage of the CCH2 herbarium consortium's own
# records (CCH2 publishes its own Darwin Core Archives to GBIF -- confirmed
# for SBBG) except UC Berkeley's Jepson Herbarium (UCJEPS), which is missing
# from GBIF entirely. That one gap was filled with a manual CCH2 portal
# export (JEPSON_DWCA_ZIP below) since CCH2 has no public query API.
#
# Usage: bash phylogeny-range.sh
# =============================================================================

#SBATCH --job-name=phylogeny_range
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=0:30:00
#SBATCH --output=phylogeny_range_%j.log

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR=/private/groups/meyerlab/esteiger/phylogeny_range
mkdir -p "$OUTDIR"

TREE_FILE="$OUTDIR/deinandra_phylogeny.nwk"
OCC_CSV="$OUTDIR/deinandra_occurrences.csv"
VILLOSA_CSV=/private/groups/meyerlab/esteiger/pheno/updatedvarsandpheno.csv
OUT_PNG="$OUTDIR/deinandra_phylogeny_range.png"
BASEMAP_DIR="$OUTDIR/basemap"
LAND_SHP="$BASEMAP_DIR/ne_10m_land.shp"

# Manual CCH2 export for UC Berkeley's Jepson Herbarium (UCJEPS) -- the one
# major CCH2-network collection not already mirrored into GBIF. Downloaded
# by hand from https://www.cch2.org/portal/ (no public query API exists);
# re-download and drop the zip in $OUTDIR if this ever needs refreshing.
JEPSON_ZIP="$OUTDIR/SymbOutput_2026-09-05_194011_DwC-A.zip"
JEPSON_DIR="$OUTDIR/jepson_dwca"
JEPSON_CSV="$JEPSON_DIR/occurrences.csv"

# GBIF genus key for Deinandra Greene (confirmed via
# https://api.gbif.org/v1/species/match?name=Deinandra&rank=GENUS)
GENUS_KEY=3106606

# Wider context window for the range map -- south-central California coast
# (Santa Barbara Co. plus adjoining San Luis Obispo / Ventura Co. margins),
# a modest zoom-out from the study bbox rather than a whole-state view.
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
# The `maps` package's built-in county coastline is too coarse right at Point
# Conception (the sharp bend the study bbox sits on) -- it cuts well inland of
# the true coast there, which makes real on-land occurrence records look like
# they're sitting offshore. Natural Earth's 10m land polygon fixes that.
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
