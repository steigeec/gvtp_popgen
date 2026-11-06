#!/bin/bash

#SBATCH --job-name=ngsF
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=100                  
#SBATCH --mem=300gb                    
#SBATCH --time=168:00:00              
#SBATCH --output=serial_test_%j.log   

cd /private/groups/meyerlab/esteiger/diversity/ngsF

# Get the input chunk of bamfiles name from command line argument
string=$1

# Check if string was provided
if [ -z "$string" ]; then
    echo "Error: Please provide a string argument"
    echo "Usage: sbatch script.sh <string>"
    exit 1
fi

GLF=${string}.glf
MAFS=${string}.mafs.gz
N_SITES=$(expr $(zcat ${MAFS} | wc -l) - 1)
N_IND=$(expr $(cat ${string} | wc -l))

pigz -d -p10 ${GLF}.gz

mkdir -p ${string}_output

run_ngsF () {
	export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH
	NUM=$1

	ngsF=/private/groups/meyerlab/esteiger/software/ngsF/ngsF
	GENOME="/private/groups/meyerlab/esteiger/genome/DEINVI_genome.fasta"

	if [ ! -f ${string}_output/500k_sites_est_${NUM}.pars ]; then
	   $ngsF --n_ind ${N_IND} --n_threads 10 --n_sites ${N_SITES} --glf ${GLF} --init_values r --approx_EM --seed $((12345 + NUM)) --min_epsilon 1e-9 --out ${string}_output/500k_sites_est_${NUM} >> ${string}_output/ngsF_log_${NUM}
	fi

	if [ ! -f ${string}_output/500k_sites_${NUM}.pars ]; then
	   $ngsF --n_ind ${N_IND} --n_threads 10 --n_sites ${N_SITES} --glf ${GLF} --init_values ${string}_output/500k_sites_est_${NUM}.pars --min_epsilon 1e-5 --out ${string}_output/500k_sites_${NUM} >> ${string}_output/ngsF_log_${NUM}
	fi
}
export -f run_ngsF
export string
export GLF
export MAFS
export N_SITES
export N_IND

parallel -j 10 run_ngsF ::: {1..100}
pigz -p10 ${GLF}

