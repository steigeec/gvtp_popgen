#!/bin/bash

#SBATCH --job-name=pop_het
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100gb
#SBATCH --time=168:00:00
#SBATCH --output=pop_het_%j.log

# SAF files are built by /private/groups/meyerlab/esteiger/gvtp_popgen/pop_SAF.sh
# This script estimates folded heterozygosity from the per-population SAF via realSFS.
# Skips populations whose .ml already exists.

SAFDIR=/private/groups/meyerlab/esteiger/SAF/population
OUTDIR=/private/groups/meyerlab/esteiger/diversity/He/population
RS=/private/groups/meyerlab/esteiger/software/angsd/misc/realSFS

POPS=(arguello conception lionhead cdfw strauss Dangermond)

for POP in "${POPS[@]}"; do
    MLFILE="${OUTDIR}/${POP}.ml"

    if [ -f "${MLFILE}" ]; then
        echo "$(date) -- ${POP}: .ml exists, skipping"
        continue
    fi

    if [ ! -f "${SAFDIR}/${POP}.saf.idx" ]; then
        echo "$(date) -- ${POP}: SAF not found in ${SAFDIR}, skipping (run pop_SAF.sh first)"
        continue
    fi

    echo "$(date) -- ${POP}: running realSFS"
    ${RS} "${SAFDIR}/${POP}.saf.idx" -fold 1 -P 15 > "${MLFILE}"
    echo "$(date) -- ${POP}: done"
done
