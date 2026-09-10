#!/bin/bash
#SBATCH --job-name=watterson
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/demography/ANGSD-Ne/logs/%x_%a.out
#SBATCH --error=/private/groups/meyerlab/esteiger/demography/ANGSD-Ne/logs/%x_%a.err
#SBATCH --array=0-5

set -euo pipefail

MISC=/private/groups/meyerlab/esteiger/software/angsd_dev/misc
SAFDIR=/private/groups/meyerlab/esteiger/SAF/population
MLDIR=/private/groups/meyerlab/esteiger/diversity/He/population
OUTDIR=/private/groups/meyerlab/esteiger/demography/ANGSD-Ne
CORES=${SLURM_CPUS_PER_TASK:-8}

POPS=(arguello cdfw conception Dangermond lionhead strauss)
POP=${POPS[$SLURM_ARRAY_TASK_ID]}

cd "$OUTDIR"

echo "[$(date)] Starting Watterson's theta pipeline for population: $POP"

# 1. Folded global SFS: reuse the one already computed for heterozygosity
#    (pop_het.sh: realSFS <pop>.saf.idx -fold 1), since it's the identical
#    SAF file + identical folded-EM call we would otherwise repeat here.
SFS="${OUTDIR}/${POP}.folded.sfs"
cp "${MLDIR}/${POP}.ml" "$SFS"
echo "[$(date)] Step 1: reused existing folded SFS from ${MLDIR}/${POP}.ml -> $SFS"

# 2. Compute per-site thetas from the SAF + folded SFS prior
echo "[$(date)] Step 2: saf2theta for $POP"
$MISC/realSFS saf2theta "${SAFDIR}/${POP}.saf.idx" -sfs "$SFS" -fold 1 -outname "${OUTDIR}/${POP}" 2> "logs/${POP}.saf2theta.log"

# 3. Genome-wide (per-contig) summary stats, incl. Watterson's theta
echo "[$(date)] Step 3: thetaStat do_stat for $POP"
$MISC/thetaStat do_stat "${OUTDIR}/${POP}.thetas.idx" 2> "logs/${POP}.thetaStat.log"

echo "[$(date)] Done: $POP"
