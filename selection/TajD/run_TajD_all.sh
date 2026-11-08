#!/bin/bash
# Tajima's D in 100-SNP windows for all five populations within one job (populations run in parallel).
#SBATCH --job-name=TajD_all
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=400gb
#SBATCH --time=168:00:00
#SBATCH --output=TajD_all_%j.log

BASE_DIR=/private/groups/meyerlab/esteiger

outdir=${BASE_DIR}/selection/TajD
sites_file=${outdir}/all_sites
regions=${outdir}/all_chrs.txt
REF=${BASE_DIR}/genome/DEINVI_genome.fasta.masked
ANGSD=${BASE_DIR}/software/angsd_dev/angsd
RS=${BASE_DIR}/software/angsd_dev/misc/realSFS
THETASTAT=${BASE_DIR}/software/angsd_dev/misc/thetaStat
TAJIMAS=${BASE_DIR}/gvtp_popgen/selection/TajD/Tajimas_SNPmidpoint.py
THREADS=4

cd ${BASE_DIR}/mapping/merged

run_tajd() {
    local name=$1
    local bam=$2
    local minind=$3
    local nindiv=$4
    local prefix="4x_${name}"

    if [[ -f "${outdir}/${prefix}.Diana.out" ]]; then
        echo "[${name}] Skipping: ${prefix}.Diana.out already exists"
        return 0
    fi

    echo "[${name}] Starting pipeline..."

    if [[ ! -f "${outdir}/${prefix}.saf.idx" ]]; then
        echo "[${name}] Running ANGSD..."
        ${ANGSD} -nThreads ${THREADS} -bam ${bam} -ref ${REF} -anc ${REF} \
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
}

bamdir=${BASE_DIR}/mapping/merged

run_tajd "strauss"    "${bamdir}/bamlist_4x_strauss.txt"    50 70 &
run_tajd "arguello"   "${bamdir}/bamlist_4x_arguello.txt"   10 14 &
run_tajd "conception" "${bamdir}/bamlist_4x_conception.txt"  11 16 &
run_tajd "lionhead"   "${bamdir}/bamlist_4x_lionhead.txt"   11 15 &
run_tajd "cdfw"       "${bamdir}/bamlist_4x_cdfw.txt"       11 16 &

wait
echo "All populations complete."
