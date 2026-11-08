#!/bin/bash
# launch_ohana_year1_strauss.sh
# Submit ohana.sh for the year1 and Strauss subsets, then (on a later run, once
# qpas has finished) pick the highest-log-likelihood seed per K into best_runs.txt.
# Usage: bash launch_ohana_year1_strauss.sh
BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/structure/OHANA/

SUBSETS=(
    "angsd_paralog_LD_year1"
    "angsd_paralog_LD_strauss"
)
# Strauss NMI plateaus after K=6 (0.628-0.652 for K=7-9, versus gains of
# 0.05-0.32 from K=2-6), so Strauss stops at K=9.
declare -A MAX_K=(
    ["angsd_paralog_LD_year1"]="10"
    ["angsd_paralog_LD_strauss"]="9"
)

# ohana.sh runs MAX_JOBS=8 qpas processes at once, each loading the full LGM, so
# --mem must cover ~8x one run's peak RSS (year1 ~9.8G, strauss ~20.7G).
# Observed peak RSS for one run: year1 ~9.8G, strauss ~20.7G.
declare -A MEM=(
    ["angsd_paralog_LD_year1"]="100G"
    ["angsd_paralog_LD_strauss"]="200G"
)

for SUBSET in "${SUBSETS[@]}"; do
    sbatch --mem=${MEM[$SUBSET]} ${BASE_DIR}/gvtp_popgen/structure/ohana.sh ${SUBSET} ${MAX_K[$SUBSET]}
done

# Select best run per-K
for SUBSET in "${SUBSETS[@]}"; do
    subdir="${SUBSET#angsd_paralog_LD_}"

    if ! ls ${subdir}/qpas_*.out 1> /dev/null 2>&1; then
        continue
    fi

    echo "Processing $subdir..."

    outfile="${subdir}/log_likelihoods_across_runs.txt"

    echo -e "K\titeration\tloglik" > "$outfile"

    for file in ${subdir}/qpas_${SUBSET}_*.out; do
        if [[ -f "$file" ]]; then
            filename=$(basename "$file")

            basename="${filename#qpas_${SUBSET}_}"
            basename="${basename%.out}"

            # e.g. "8_3" (K_seed)
            K=$(echo "$basename" | cut -d'_' -f1)
            iteration=$(echo "$basename" | cut -d'_' -f2)

            loglik=$(grep "^500" "$file" | cut -f3)

            if [[ -z "$loglik" ]]; then
                echo "Warning: No line starting with 500 found in $file, skipping" >&2
                continue
            fi

            echo -e "${K}\t${iteration}\t${loglik}" >> "$outfile"
        fi
    done

    # Best (max log-likelihood) seed per K
    if [[ -f "$outfile" ]] && [[ $(wc -l < "$outfile") -gt 1 ]]; then
        echo "Best runs for $subdir:"
        awk 'NR==1 {next}
        {
            if (!($1 in best) || $3 > best[$1]) {
                best[$1] = $3
                line[$1] = $0
            }
        }
        END {
            n = asorti(line, sorted)
            for (i=1; i<=n; i++) print line[sorted[i]]
        }' "$outfile" > "${subdir}/best_runs.txt"

        cat "${subdir}/best_runs.txt"
        echo ""
    fi
done
