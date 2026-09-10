#!/bin/bash
# Submit one SLURM job per population with memory scaled to population size.

SCRIPT=/private/groups/meyerlab/esteiger/gvtp_popgen/selection/TajD/run_TajD.sh
LOGDIR=/private/groups/meyerlab/esteiger/gvtp_popgen/selection/TajD

#                name        minInd  nIndiv  threads  mem
populations=(
    "strauss      50       70       15       400gb"
    "arguello     10       14        4        50gb"
    "conception   11       16        4        50gb"
    "lionhead     11       15        4        50gb"
    "cdfw         11       16        4        50gb"
)

for pop in "${populations[@]}"; do
    read -r name minind nindiv threads mem <<< "$pop"
    echo "Submitting ${name} (minInd=${minind}, nIndiv=${nindiv}, threads=${threads}, mem=${mem})..."
    sbatch --job-name="TajD_${name}" \
           --ntasks="${threads}" \
           --mem="${mem}" \
           --output="${LOGDIR}/TajD_${name}_%j.log" \
           "${SCRIPT}" "${name}" "${minind}" "${nindiv}"
done
