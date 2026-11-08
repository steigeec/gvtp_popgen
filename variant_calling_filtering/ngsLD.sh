#!/bin/bash
# ngsLD on one per-scaffold GL file (>4x Strauss individuals).
# Usage: ngsLD.sh <scaffold.tsv.gz>  (run by ngsLD_launch.sh)

BASE_DIR=/private/groups/meyerlab/esteiger

TSV=$1
prefix="./"
TSV=${TSV#"$prefix"}
suffix=".tsv.gz"
string=${TSV%"$suffix"}

cd ${BASE_DIR}/variant_calling/ngsLD

NGSLD="${BASE_DIR}/software/ngsLD/ngsLD"
BAMLIST="${BASE_DIR}/mapping/merged/bamlist_4x_strauss.txt"
SITES_DB="${BASE_DIR}/variant_calling/angsd/postParalog_MAF.sites"

N_INDIVIDUALS=$(wc -l < "$BAMLIST")
N_SITES=$(zcat "$TSV" | wc -l)
SITES="${string}.sites"

if [ -f "${string}_prune.ld" ]; then
    echo "Output already exists, skipping ${string}"
    exit 0
fi

# positions for this scaffold, sorted by chr then numeric position
grep "${string}" "$SITES_DB" | sort -k1,1 -k2,2n > "$SITES"

# LD for all site pairs within 300 kb
# (300 kb chosen because some apparent structural variation spans this distance)
$NGSLD \
    --geno "$TSV" \
    --probs \
    --n_ind "$N_INDIVIDUALS" \
    --n_sites "$N_SITES" \
    --pos "$SITES" \
    --max_kb_dist 300 \
    --max_snp_dist 0 \
    --n_threads 30 \
    --out "${string}_prune.ld"
