#!/bin/bash
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --time=168:00:00

# Tajima's D in 100-SNP windows for one population (ANGSD, realSFS, thetaStat, Tajimas_SNPmidpoint.py).
# Usage: sbatch --job-name=TajD_<name> --ntasks=<threads> --mem=<mem>
#        --output=TajD_<name>_%j.log run_TajD.sh <name> <minInd> <nIndiv>

BASE_DIR=/private/groups/meyerlab/esteiger

name=$1
minind=$2
nindiv=$3
prefix="4x_${name}"

outdir=${BASE_DIR}/selection/TajD
sites_file=${outdir}/all_sites
regions=${outdir}/all_chrs.txt
REF=${BASE_DIR}/genome/DEINVI_genome.fasta.masked
ANGSD=${BASE_DIR}/software/angsd_dev/angsd
RS=${BASE_DIR}/software/angsd_dev/misc/realSFS
THETASTAT=${BASE_DIR}/software/angsd_dev/misc/thetaStat
TAJIMAS=${BASE_DIR}/gvtp_popgen/selection/TajD/Tajimas_SNPmidpoint.py
BAM=${BASE_DIR}/mapping/merged/bamlist_4x_${name}.txt
THREADS=${SLURM_NTASKS:-4}

cd ${BASE_DIR}/mapping/merged

if [[ -f "${outdir}/${prefix}.Diana.out" ]]; then
    echo "[${name}] Skipping: ${prefix}.Diana.out already exists"
    exit 0
fi

echo "[${name}] Starting pipeline (threads=${THREADS}, minInd=${minind}, nIndiv=${nindiv})..."

if [[ ! -f "${outdir}/${prefix}.saf.idx" ]]; then
    echo "[${name}] Running ANGSD..."
    ${ANGSD} -nThreads ${THREADS} -bam ${BAM} -ref ${REF} -anc ${REF} \
        -out ${outdir}/${prefix} \
        -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
        -minMapQ 20 -minQ 20 -minInd ${minind} -setMinDepthInd 1 \
        -setMaxDepthInd 25 -doSaf 1 -rf ${regions} -sites ${sites_file} -GL 2
else
    echo "[${name}] Skipping ANGSD: ${prefix}.saf.idx already exists"
fi

if [[ ! -f "${outdir}/${prefix}.ml" ]]; then
    echo "[${name}] Running realSFS..."
    ${RS} ${outdir}/${prefix}.saf.idx -fold 1 -P ${THREADS} \
        > ${outdir}/${prefix}.ml
else
    echo "[${name}] Skipping realSFS: ${prefix}.ml already exists"
fi

if [[ ! -f "${outdir}/${prefix}.thetas.idx" ]]; then
    echo "[${name}] Running saf2theta..."
    ${RS} saf2theta ${outdir}/${prefix}.saf.idx \
        -sfs ${outdir}/${prefix}.ml -outname ${outdir}/${prefix}
else
    echo "[${name}] Skipping saf2theta: ${prefix}.thetas.idx already exists"
fi

if [[ ! -f "${outdir}/${prefix}.thetas.persite.txt" ]]; then
    echo "[${name}] Running thetaStat..."
    ${THETASTAT} print ${outdir}/${prefix}.thetas.idx \
        > ${outdir}/${prefix}.thetas.persite.txt
else
    echo "[${name}] Skipping thetaStat: ${prefix}.thetas.persite.txt already exists"
fi

echo "[${name}] Running Tajimas_SNPmidpoint.py..."
${TAJIMAS} ${outdir}/${prefix}.thetas.persite.txt \
    ${outdir}/${prefix}.Diana.out 100 ${nindiv}

echo "[${name}] Done."
