#!/bin/bash
# Sensitivity check for population pairwise Fst: recompute all 15 pairs with
# every population downsampled to n = 4 (Dangermond's sample size), 3 replicates.
# Dangermond keeps its existing n = 4 SAF. Fst is computed as in
# structure/calculate_Fst_pairwisepop.sh.
#   random    4 individuals drawn at random from the population's >= 4x bamlist
#   covmatch  4 individuals each within 1.5x of one Dangermond coverage (nearest if none)
#
# Pipeline (SLURM dependency chain):
#   fst_downsample_saf.sh   array, one task per population x replicate
#   fst_downsample_pair.sh  array, one task per pair x replicate
#   collect                 one table: <OUTDIR>/pairwise_Fst_downsampled.txt
# Outputs: wind/fst_downsampled (random), wind/fst_downsampled_covmatch.
# Environment: python3, SLURM; angsd/realSFS paths are set in the task scripts.
#
# Usage: bash fst_downsample_launch.sh random|covmatch

set -euo pipefail

export BASE_DIR=/private/groups/meyerlab/esteiger

export BAMDIR=${BASE_DIR}/mapping/merged
DESIGN="${1:?usage: fst_downsample_launch.sh random|covmatch}"
case "${DESIGN}" in
  random)   export OUTDIR=${BASE_DIR}/wind/fst_downsampled ;;
  covmatch) export OUTDIR=${BASE_DIR}/wind/fst_downsampled_covmatch ;;
  *) echo "unknown design ${DESIGN}"; exit 1 ;;
esac
export FULL_SAFDIR=${BASE_DIR}/SAF/population
SCRIPTS_DIR=${BASE_DIR}/gvtp_popgen/wind
N_DRAW=4
REPS=(1 2 3)
mkdir -p "${OUTDIR}/bamlists" "${OUTDIR}/saf" "${OUTDIR}/pairs" "${OUTDIR}/logs"

declare -A FULL_LISTS=(
  [arguello]=bamlist_4x_arguello.txt
  [conception]=bamlist_4x_conception.txt
  [lionhead]=bamlist_4x_lionhead.txt
  [cdfw]=bamlist_4x_cdfw.txt
  [strauss]=bamlist_4x_strauss_top16.txt
)
POPS=(conception cdfw arguello strauss lionhead Dangermond)

declare -A CLUSTERS=(
  [arguello]=VAFB_Arguello [conception]=Point_Conception [lionhead]=VAFB_Lion_Head
  [cdfw]=CDFW [strauss]=Strauss
)

# ---------- n = 4 draws (RNG seeded with "<pop>-rep<r>") ----------
: > "${OUTDIR}/saf_tasks.txt"
for rep in "${REPS[@]}"; do
  for pop in "${!FULL_LISTS[@]}"; do
    list="${OUTDIR}/bamlists/${pop}_rep${rep}.txt"
    if [ "${DESIGN}" = random ]; then
      python3 -c "import random, sys; ids = [l.strip() for l in open(sys.argv[1]) if l.strip()]; \
print('\n'.join(sorted(random.Random(sys.argv[2]).sample(ids, int(sys.argv[3])))))" \
        "${BAMDIR}/${FULL_LISTS[$pop]}" "${pop}-rep${rep}" "${N_DRAW}" | sed "s|^|${BAMDIR}/|" > "${list}"
    else
      python3 - "${CLUSTERS[$pop]}" "${pop}-rep${rep}" "${BAMDIR}" > "${list}" <<'PY'
import csv, math, os, random, sys
cluster, seed, bamdir = sys.argv[1:4]
rng = random.Random(seed)
# 13 individuals with untrustworthy coordinates
bad = {"4548", "4549", "4550", "4551", "4552", "4553", "4554", "4555", "4556", "4557", "4559", "4561", "4562"}
rows = list(csv.DictReader(open(os.path.join(os.environ["BASE_DIR"], "field_data.txt")), delimiter="\t"))
cov = {r["SampleID"]: float(r["Cov"]) for r in rows if r["Cov"] not in ("", "NA")}
dang = [l.strip().replace("_rmdup.bam", "") for l in open(f"{bamdir}/bamlist_Dangermond.txt") if l.strip()]
pool = [r["SampleID"] for r in rows if r["Clusters"] == cluster and r["SampleID"] not in bad
        and r["SampleID"] in cov and os.path.exists(f"{bamdir}/{r['SampleID']}_rmdup.bam")]
chosen = []
for target in sorted(cov[d] for d in dang):
    left = [i for i in pool if i not in chosen]
    ok = [i for i in left if abs(math.log(cov[i] / target)) <= math.log(1.5)]
    pick = rng.choice(sorted(ok)) if ok else min(left, key=lambda i: abs(math.log(cov[i] / target)))
    if not ok:
        print(f"{cluster} {seed}: no individual within 1.5x of {target}x; nearest {pick} ({cov[pick]}x)", file=sys.stderr)
    chosen.append(pick)
print("\n".join(f"{bamdir}/{i}_rmdup.bam" for i in chosen))
PY
    fi
    echo "${pop} ${rep}" >> "${OUTDIR}/saf_tasks.txt"
  done
done

# ---------- All pairs x replicates (pair order as in calculate_Fst_pairwisepop.sh) ----------
: > "${OUTDIR}/pair_tasks.txt"
for rep in "${REPS[@]}"; do
  for (( i=0; i<${#POPS[@]}-1; i++ )); do
    for (( j=i+1; j<${#POPS[@]}; j++ )); do
      echo "${POPS[$i]} ${POPS[$j]} ${rep}" >> "${OUTDIR}/pair_tasks.txt"
    done
  done
done
N_SAF=$(wc -l < "${OUTDIR}/saf_tasks.txt"); N_PAIR=$(wc -l < "${OUTDIR}/pair_tasks.txt")
echo "$(date) -- ${N_SAF} downsampled SAFs, ${N_PAIR} pairwise comparisons"

SAF_JID=$(sbatch --array=1-${N_SAF} --export=ALL --output="${OUTDIR}/logs/saf_%A_%a.log" \
  "${SCRIPTS_DIR}/fst_downsample_saf.sh" | awk '{print $NF}')
PAIR_JID=$(sbatch --array=1-${N_PAIR} --export=ALL --dependency=afterok:${SAF_JID} \
  --output="${OUTDIR}/logs/pair_%A_%a.log" "${SCRIPTS_DIR}/fst_downsample_pair.sh" | awk '{print $NF}')
COL_JID=$(sbatch --job-name=fstds_collect --partition=short --ntasks=1 --mem=2gb --time=0:30:00 \
  --dependency=afterok:${PAIR_JID} \
  --output="${OUTDIR}/logs/collect_%j.log" \
  --wrap="{ printf 'pop1\tpop2\treplicate\tFst_unweighted\tFst_weighted\n'; cat ${OUTDIR}/pairs/*.stat; } > ${OUTDIR}/pairwise_Fst_downsampled.txt" \
  | awk '{print $NF}')
echo "Pipeline: SAF ${SAF_JID} -> pairs ${PAIR_JID} -> collect ${COL_JID}"
