#!/bin/bash --login

#SBATCH --job-name=FEEMS
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=1                  
#SBATCH --mem=300gb                     # 5gb per job * 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=/private/groups/meyerlab/esteiger/migration/FEEMS/FEEMS_year1.log   # Standard output and error log

# Conda env: feems_e
#   conda create -n feems_e -c conda-forge -c bioconda python=3.8.10 \
#     numpy=1.22.3 scipy=1.5.0 scikit-learn=0.23.1 matplotlib=3.2.2 \
#     networkx=2.4 pyproj=2.6.1 shapely=1.7.1 fiona cartopy=0.18.0 \
#     suitesparse=5.7.2 scikit-sparse=0.4.4 pandas-plink msprime=1.0.0 \
#     statsmodels=0.14.0
#   pip install git+https://github.com/NovembreLab/feems@c08a179   # FEEMS 1.0.1
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