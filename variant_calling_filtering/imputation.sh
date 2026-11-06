#!/bin/bash

#SBATCH --job-name=imputation
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=300gb
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=serial_test_%j.log   # Standard output and error log

export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH

VCF1="$1"
if [[ -z "$VCF1" ]]; then
    echo "Usage: $0 <input.bcf>" >&2
    exit 1
fi

BASENAME=$(basename "$VCF1" .bcf)
OUT_PREFIX="${BASENAME}_GLimp"

BGL=/private/groups/meyerlab/esteiger/software/beagle.r1399.jar

cd "$(dirname "$VCF1")"

#  imputation (BEAGLE r1399 cannot read BCF directly; pipe through bcftools)
bcftools view "$VCF1" | java -Xmx15000m -jar $BGL nthreads=50 gl=/dev/stdin out="$OUT_PREFIX" impute=true

# Filter by imputation quality (DR2 > 0.8)
bcftools view -i 'INFO/DR2>0.8' "${OUT_PREFIX}.vcf.gz" > "${OUT_PREFIX}.vcf"
