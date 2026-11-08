#!/bin/bash
set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

# Union five GO annotation sources into one gene_id/go_term table:
#   eggNOG-mapper, Viridiplantae scope (convert_eggnog_to_geneid_go.py)
#   BLAST-hit harvest, no_go and has_go modes (harvest_blast_go.py)
#   1:1 Arabidopsis orthologs (build_tair_ortholog_go.py)
#   1:1 Helianthus annuus orthologs (build_ortholog_go.py)
# Lettuce and tomato are not used: their BioMart GO annotation is IEA only.
#
# Usage: bash merge_go_sources.sh

GENOME_DIR=${BASE_DIR}/genome
EGGNOG="${GENOME_DIR}/annotated_sequences_eggnog_viridiplantae_geneid_go.tsv"
HARVEST="${GENOME_DIR}/annotated_sequences_blastharvest_geneid_go.tsv"
HARVEST_SUPPLEMENT="${GENOME_DIR}/annotated_sequences_blastharvest_supplement_geneid_go.tsv"
ORTHOFINDER_TAIR="${GENOME_DIR}/annotated_sequences_orthofinder_tair_geneid_go.tsv"
ORTHOFINDER_HANNUUS="${GENOME_DIR}/annotated_sequences_orthofinder_hannuus_geneid_go.tsv"
COMBINED="${GENOME_DIR}/annotated_sequences_combined_geneid_go.tsv"

for f in "$EGGNOG" "$HARVEST" "$HARVEST_SUPPLEMENT" "$ORTHOFINDER_TAIR" "$ORTHOFINDER_HANNUUS"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: missing or empty input file: $f" >&2
        exit 1
    fi
done

echo -e "gene_id\tgo_term" > "$COMBINED"
tail -n +2 -q "$EGGNOG" "$HARVEST" "$HARVEST_SUPPLEMENT" "$ORTHOFINDER_TAIR" "$ORTHOFINDER_HANNUUS" | sort -u >> "$COMBINED"

echo "Wrote $(tail -n +2 "$COMBINED" | wc -l) gene-GO pairs across " \
     "$(tail -n +2 "$COMBINED" | cut -f1 | sort -u | wc -l) genes to $COMBINED"
echo "  (eggNOG alone: $(tail -n +2 "$EGGNOG" | cut -f1 | sort -u | wc -l) genes, " \
     "no_go harvest alone: $(tail -n +2 "$HARVEST" | cut -f1 | sort -u | wc -l) genes, " \
     "has_go harvest supplement alone: $(tail -n +2 "$HARVEST_SUPPLEMENT" | cut -f1 | sort -u | wc -l) genes, " \
     "OrthoFinder/TAIR ortholog alone: $(tail -n +2 "$ORTHOFINDER_TAIR" | cut -f1 | sort -u | wc -l) genes, " \
     "OrthoFinder/sunflower ortholog alone: $(tail -n +2 "$ORTHOFINDER_HANNUUS" | cut -f1 | sort -u | wc -l) genes)"
