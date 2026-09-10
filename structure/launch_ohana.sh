#!/bin/bash

cd /private/groups/meyerlab/esteiger/structure/OHANA/

SUBSETS=(
    #"angsd_paralog_LD_year1"
    #"angsd_paralog_LD_strauss"
    "angsd_paralog_LD_CDFW"
    "angsd_paralog_LD_other"
    "angsd_paralog_LD_VAFB_Lions_Head"
)

for SUBSET in "${SUBSETS[@]}"; do
    sbatch /private/groups/meyerlab/esteiger/gvtp_popgen/structure/ohana.sh ${SUBSET}
done

# Select best run per-K
for subdir in */; do
    # Remove trailing slash
    subdir=${subdir%/}

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
    for file in ${subdir}/qpas_angsd_paralog_LD_${subdir}_*.out; do
        if [[ -f "$file" ]]; then
            # Extract just the filename
            filename=$(basename "$file")

            # Remove the prefix "qpas_angsd_paralog_LD_${subdir}_" and suffix ".out"
            basename="${filename#qpas_angsd_paralog_LD_${subdir}_}"
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
