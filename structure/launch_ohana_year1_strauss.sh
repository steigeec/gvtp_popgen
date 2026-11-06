#!/bin/bash

cd /private/groups/meyerlab/esteiger/structure/OHANA/

SUBSETS=(
    "angsd_paralog_LD_year1"
    "angsd_paralog_LD_strauss"
)
# Strauss's NMI plateaus after K=6 (0.628-0.652 across K=7-9, an order of
# magnitude smaller than the 0.05-0.32 gains seen from K=2-6), so K=10 was
# dropped for strauss as not worth the compute. year1's K exploration is
# unaffected.
declare -A MAX_K=(
    ["angsd_paralog_LD_year1"]="10"
    ["angsd_paralog_LD_strauss"]="9"
)

# ohana.sh runs MAX_JOBS=8 qpas processes concurrently, and each one loads the
# full LGM independently, so --mem needs to cover ~8x a single run's peak RSS.
# Observed peak RSS for one run: year1 ~9.8G, strauss ~20.7G.
declare -A MEM=(
    ["angsd_paralog_LD_year1"]="100G"
    ["angsd_paralog_LD_strauss"]="200G"
)

for SUBSET in "${SUBSETS[@]}"; do
    sbatch --mem=${MEM[$SUBSET]} /private/groups/meyerlab/esteiger/gvtp_popgen/structure/ohana.sh ${SUBSET} ${MAX_K[$SUBSET]}
done

# Select best run per-K
for SUBSET in "${SUBSETS[@]}"; do
    subdir="${SUBSET#angsd_paralog_LD_}"

    # Check if this directory has qpas files
    if ! ls ${subdir}/qpas_*.out 1> /dev/null 2>&1; then
        continue
    fi

    echo "Processing $subdir..."

    # Output file for this subdirectory
    outfile="${subdir}/log_likelihoods_across_runs.txt"

    # Clear/create output file with header
    echo -e "K\titeration\tloglik" > "$outfile"

    # Process all qpas .out files in this subdirectory
    for file in ${subdir}/qpas_${SUBSET}_*.out; do
        if [[ -f "$file" ]]; then
            # Extract just the filename
            filename=$(basename "$file")

            # Remove the prefix "qpas_${SUBSET}_" and suffix ".out"
            basename="${filename#qpas_${SUBSET}_}"
            basename="${basename%.out}"

            # Now basename should be like "8_3"
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

    # Choose the max log likelihood iteration per K
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
