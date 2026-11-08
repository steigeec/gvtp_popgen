#!/bin/bash
# ohana.sh
# OHANA qpas for K=2..MAX_K with 10 seeds each, MAX_JOBS runs in parallel.
# Usage: sbatch --mem=<mem> ohana.sh <subset> [max_K]   (see launch_ohana_year1_strauss.sh)
#SBATCH --job-name=ohana
#SBATCH --output=ohana_strauss.out
#SBATCH --time=168:00:00            # Wall time (adjust as needed)
#SBATCH --cpus-per-task=8           # Must stay in sync with MAX_JOBS below
#SBATCH --nodes=1
# --mem: qpas loads the whole LGM per process and MAX_JOBS run at once; pass --mem per subset
#SBATCH --mem=32G
#SBATCH --partition=long       # Partition (adjust as needed)

BASE_DIR=/private/groups/meyerlab/esteiger

SUBSET=$1
MAX_K=${2:-8}
MAX_JOBS=8

export PYTHONUNBUFFERED=1
set -x

CON="${BASE_DIR}/software/ohana/bin/convert"
BEAGLE_DIR="${BASE_DIR}/variant_calling/angsd/paralog_filtered"
QPAS="${BASE_DIR}/software/ohana/bin/qpas"
NEM="${BASE_DIR}/software/ohana/bin/nemeco"
LGM=${SUBSET}.lgm
BEAGLE="${BEAGLE_DIR}/${SUBSET}.beagle"

if [ ! -f ${LGM} ]; then
  ${CON} bgl2lgm ${BEAGLE} ${LGM}
fi

cd ${BASE_DIR}/structure/OHANA

SUBDIR="${SUBSET#angsd_paralog_LD_}"
mkdir -p ${SUBDIR}

for i in $(seq 2 ${MAX_K}); do
 for n in {1..10}; do

    QO="${SUBDIR}/q_${SUBSET}_${i}_${n}.matrix"
    FO="${SUBDIR}/f_${SUBSET}_${i}_${n}.matrix"
    CO="${SUBDIR}/c_${SUBSET}_${i}_${n}.matrix"
    OUT="${SUBDIR}/qpas_${SUBSET}_${i}_${n}.out"

    if [ -f ${OUT} ]; then
        continue
    fi

    if [ ! -f ${QO} ]; then
        echo "Running qpas for $SUBSET, K=$i, seed=$n"
        $QPAS $LGM -k ${i} -qo $QO -fo $FO -mi 500 -s ${n} > ${OUT} &
    fi

    # Cap concurrent qpas processes at MAX_JOBS
    while [ "$(jobs -r -p | wc -l)" -ge ${MAX_JOBS} ]; do
        wait -n
    done
 done
done

# Wait for background qpas runs
wait
