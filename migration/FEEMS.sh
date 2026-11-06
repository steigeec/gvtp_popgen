#!/bin/bash --login

#SBATCH --job-name=FEEMS
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=1                  
#SBATCH --mem=300gb                     # 5gb per job * 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=/private/groups/meyerlab/esteiger/migration/FEEMS/FEEMS_year1.log   # Standard output and error log

BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/migration/FEEMS

# ---- Per-year plink inputs (year{1,2,3}_strauss_maf) ----
# LD-pruned VCF (see variant_calling_filtering/make_LD_sites.sh) subset to each
# year's samples; --mac 1 drops sites monomorphic within the year.
LD_VCF="${BASE_DIR}/variant_calling/ngsLD/angsd_paralog_LD_GLimp.vcf"
YEAR_LIST_DIR="${BASE_DIR}/mapping/merged"
FIELD_DATA="${BASE_DIR}/field_data.txt"
PLINK="${BASE_DIR}/software/plink"

for N in 1 2 3; do
    YEAR="year${N}_strauss_maf"
    bcftools view -S "${YEAR_LIST_DIR}/strauss_year${N}.txt" "$LD_VCF" -Oz -o "${YEAR}.vcf.gz"
    "$PLINK" --vcf "${YEAR}.vcf.gz" --allow-extra-chr --double-id --mac 1 --make-bed --out "$YEAR"

    # coord file: Long, Lat per sample, same row order as ${YEAR}.fam
    awk -v OFS='\t' '
        NR==FNR { if (FNR>1) { lat[$1]=$4; lon[$1]=$5 }; next }
        {
            id=$2; sub(/_rmdup\.bam$/,"",id)
            if (!(id in lat)) { print "ERROR: no field_data.txt entry for sample " id > "/dev/stderr"; err=1; next }
            print lon[id], lat[id]
        }
        END { if (err) exit 1 }
    ' "$FIELD_DATA" "${YEAR}.fam" > "${YEAR}.coord"
done

source /private/home/esteiger/miniconda3/etc/profile.d/conda.sh
conda activate feems_e || { echo "Failed to activate Conda environment"; exit 1; }
echo "conda path is $CONDA_PREFIX"

YEARS=(
"year1_strauss_maf"
"year2_strauss_maf"
"year3_strauss_maf"
)

for YEAR in "${YEARS[@]}"; do
    $CONDA_PREFIX/bin/python ${BASE_DIR}/gvtp_popgen/migration/FEEMS.py ${YEAR}
done