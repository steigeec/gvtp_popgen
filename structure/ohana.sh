#!/bin/bash

#SBATCH --job-name=ohana    
#SBATCH --output=ohana_strauss.out  
#SBATCH --time=168:00:00            # Wall time (adjust as needed)
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1    	
#SBATCH --mem=100G                  # Memory per node (adjust as needed)
#SBATCH --partition=long       # Partition (adjust as needed)

SUBSET=$1

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

for i in {2..8}; do
 for n in {1..10}; do

    QO="${SUBDIR}/q_${SUBSET}_${i}_${n}.matrix"
    FO="${SUBDIR}/f_${SUBSET}_${i}_${n}.matrix"
    CO="${SUBDIR}/c_${SUBSET}_${i}_${n}.matrix"
    OUT="${SUBDIR}/qpas_${SUBSET}_${i}_${n}.out"
  
    if [ -f ${OUT} ]; then
        continue
    fi

    if [ ! -f ${QO} ]; then   
        echo "Running qpas for $SUBSET and K=$i"
        $QPAS $LGM -k ${i} -qo $QO -fo $FO -mi 500 -s ${n} > ${OUT} 
    fi

#    if [ ! -f ${CO} ]; then   
#        echo "Estimating component covariances (population covariances) for $SUBSET and K=$i"
#        ${NEM} $LGM $FO -co $CO -mi 500
#    fi
 done
done
