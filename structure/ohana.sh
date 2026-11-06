#!/bin/bash

#SBATCH --job-name=ohana
#SBATCH --output=ohana_strauss.out
#SBATCH --time=168:00:00            # Wall time (adjust as needed)
#SBATCH --cpus-per-task=8           # Must stay in sync with MAX_JOBS below
#SBATCH --nodes=1
#SBATCH --mem=32G                   # Default; callers should pass --mem on the
                                     # sbatch command line to override per subset
                                     # (qpas loads the whole LGM per process, and
                                     # MAX_JOBS run concurrently, so the right
                                     # value scales with both subset size and
                                     # MAX_JOBS)
#SBATCH --partition=long       # Partition (adjust as needed)

SUBSET=$1
MAX_K=${2:-8}
MAX_JOBS=8

export PYTHONUNBUFFERED=1
set -x

CON="/private/groups/meyerlab/esteiger/software/ohana/bin/convert"
BEAGLE_DIR="/private/groups/meyerlab/esteiger/variant_calling/angsd/paralog_filtered"
QPAS="/private/groups/meyerlab/esteiger/software/ohana/bin/qpas"
NEM="/private/groups/meyerlab/esteiger/software/ohana/bin/nemeco"
LGM=${SUBSET}.lgm
BEAGLE="${BEAGLE_DIR}/${SUBSET}.beagle"

if [ ! -f ${LGM} ]; then
  ${CON} bgl2lgm ${BEAGLE} ${LGM}
fi

cd /private/groups/meyerlab/esteiger/structure/OHANA

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

#    if [ ! -f ${CO} ]; then
#        echo "Estimating component covariances (population covariances) for $SUBSET and K=$i"
#        ${NEM} $LGM $FO -co $CO -mi 500
#    fi

    # Cap concurrently-running qpas processes at MAX_JOBS
    while [ "$(jobs -r -p | wc -l)" -ge ${MAX_JOBS} ]; do
        wait -n
    done
 done
done

# Don't let the job exit while background qpas runs are still in flight
wait
