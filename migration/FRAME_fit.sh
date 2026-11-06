#!/bin/bash --login

#SBATCH --job-name=FRAME_fit
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300gb
#SBATCH --time=168:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/migration/FEEMS/FRAME_pooled.log

# FRAME on all Tranquillon Mtn individuals, all years pooled (strauss_all_maf,
# 736 individuals), on the FEEMS grid. lambda_m is cross-validated; append a
# value to the python call to fix it.
# The pooled plink input is built from the LD-pruned VCF as in FEEMS.sh,
# once; skipped if it already exists.
# Env: frame_e conda env (python 3.11.9, numpy 1.26.4, scipy 1.11.4, MKL BLAS).
# PYTHONNOUSERSITE=1 keeps packages in ~/.local out of the environment.

set -euo pipefail

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/migration/FEEMS

PREFIX=strauss_all_maf
LD_VCF="${BASE_DIR}/variant_calling/ngsLD/angsd_paralog_LD_GLimp.vcf"
YEAR_LIST_DIR="${BASE_DIR}/mapping/merged"
FIELD_DATA="${BASE_DIR}/field_data.txt"
PLINK="${BASE_DIR}/software/plink"

if [ ! -s "${PREFIX}.bed" ] || [ ! -s "${PREFIX}.coord" ]; then
    cat "${YEAR_LIST_DIR}"/strauss_year{1,2,3}.txt > "${PREFIX}_samples.txt"
    bcftools view --threads 8 -S "${PREFIX}_samples.txt" "$LD_VCF" -Oz -o "${PREFIX}.vcf.gz"
    "$PLINK" --vcf "${PREFIX}.vcf.gz" --allow-extra-chr --double-id --mac 1 --make-bed --out "$PREFIX"

    # coord file: Long, Lat per sample, same row order as ${PREFIX}.fam
    awk -v OFS='\t' '
        NR==FNR { if (FNR>1) { lat[$1]=$4; lon[$1]=$5 }; next }
        {
            id=$2; sub(/_rmdup\.bam$/,"",id)
            if (!(id in lat)) { print "ERROR: no field_data.txt entry for sample " id > "/dev/stderr"; err=1; next }
            print lon[id], lat[id]
        }
        END { if (err) exit 1 }
    ' "$FIELD_DATA" "${PREFIX}.fam" > "${PREFIX}.coord"
fi
echo "${PREFIX}: $(wc -l < "${PREFIX}.fam") individuals, $(wc -l < "${PREFIX}.bim") sites"

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
set +u; conda activate frame_e || { echo "Failed to activate Conda environment"; exit 1; }; set -u
echo "conda path is $CONDA_PREFIX"
export PYTHONNOUSERSITE=1

"$CONDA_PREFIX/bin/python" ${BASE_DIR}/gvtp_popgen/migration/FRAME_fit.py "$PREFIX"
