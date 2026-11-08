#!/bin/bash
#SBATCH --job-name=watterson_summary
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/private/groups/meyerlab/esteiger/demography/ANGSD-Ne/logs/%x.out
#SBATCH --error=/private/groups/meyerlab/esteiger/demography/ANGSD-Ne/logs/%x.err

# Sum Watterson's theta and sites over contigs per population and convert to
# Ne = theta_W / (4 mu). Writes watterson_theta_Ne_summary.tsv.

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

OUTDIR=${BASE_DIR}/demography/ANGSD-Ne
MU=1e-8   # per site per generation; Asteraceae rate (Strasburg & Rieseberg 2008,
          # Helianthus annuus/petiolaris; Senecio chrysanthemifolius, ~1e-8)

POPS=(arguello cdfw conception Dangermond lionhead strauss)
OUT="${OUTDIR}/watterson_theta_Ne_summary.tsv"

echo -e "population\tsum_tW\tsum_nSites\ttheta_W_per_site\tNe" > "$OUT"

cd "$OUTDIR"

for POP in "${POPS[@]}"; do
    PESTPG="${POP}.thetas.idx.pestPG"
    if [[ ! -f "$PESTPG" ]]; then
        echo "WARNING: $PESTPG not found, skipping $POP" >&2
        continue
    fi
    # pestPG columns (no header token on col1): ...  tW  tP  tF  tH  tL  Tajima  fuf  fud  fayh  zeng  nSites
    # tW = per-window sum of Watterson's theta; nSites = sites with data in window.
    # windows are whole contigs (do_stat run without -win/-step); sum for genome-wide totals
    read sum_tW sum_nSites <<< "$(awk 'NR>1 {tw+=$4; ns+=$14} END{print tw, ns}' "$PESTPG")"

    theta_w_per_site=$(awk -v tw="$sum_tW" -v ns="$sum_nSites" 'BEGIN{printf "%.10g", tw/ns}')
    ne=$(awk -v tw="$theta_w_per_site" -v mu="$MU" 'BEGIN{printf "%.2f", tw/(4*mu)}')

    echo -e "${POP}\t${sum_tW}\t${sum_nSites}\t${theta_w_per_site}\t${ne}" >> "$OUT"
done

echo "Wrote summary to $OUT"
cat "$OUT"
