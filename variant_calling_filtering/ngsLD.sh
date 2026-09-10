#!/bin/bash

TSV=$1
prefix="./"
TSV=${TSV#"$prefix"}
suffix=".tsv.gz"
string=${TSV%"$suffix"}

cd /private/groups/meyerlab/esteiger/variant_calling/ngsLD

NGSLD="/private/groups/meyerlab/esteiger/software/ngsLD/ngsLD"
BAMLIST="/private/groups/meyerlab/esteiger/mapping/merged/bamlist_4x_strauss.txt"
SITES_DB="/private/groups/meyerlab/esteiger/variant_calling/angsd/postParalog_MAF.sites"

N_INDIVIDUALS=$(wc -l < "$BAMLIST")
N_SITES=$(zcat "$TSV" | wc -l)
SITES="${string}.sites"

if [ -f "${string}_prune.ld" ]; then
    echo "Output already exists, skipping ${string}"
    exit 0
fi

# Create sorted (lexicographic chr, numeric pos) sites file for this chunk
grep "${string}" "$SITES_DB" | sort -k1,1 -k2,2n > "$SITES"

# Calculate LD for all sites within 300kb
# (300kb chosen because some apparent structural variation spans this distance)
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
