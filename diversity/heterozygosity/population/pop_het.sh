#!/bin/bash

#SBATCH --job-name=pop_het
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100gb
#SBATCH --time=168:00:00
#SBATCH --output=pop_het_%j.log

# pop_het.sh
# Folded SFS per population from the SAF files built by pop_SAF.sh (realSFS).
# Populations with an existing .ml are skipped.

BASE_DIR=/private/groups/meyerlab/esteiger

SAFDIR=${BASE_DIR}/SAF/population
OUTDIR=${BASE_DIR}/diversity/He/population
RS=${BASE_DIR}/software/angsd/misc/realSFS

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
