#!/bin/bash
#SBATCH -p medium -t 8:00:00 --mem 32G -c 16 -J cenh3_dnds
#SBATCH -o /private/groups/meyerlab/esteiger/inversions/centromere/cenh3_dnds/cenh3_dnds.%j.log
# Run cenh3_dnds.py with 16 processes. Usage: sbatch run_cenh3_dnds.sh
# Conda env: dnds (BASE_DIR/software/miniconda3/envs/dnds)
#   conda create -n dnds -c conda-forge -c bioconda paml=4.10.10 mafft=7.526 \
#     pal2nal=14.1 miniprot=0.18
#   (cenh3_dnds.py calls the tools by path; no activation needed)
BASE_DIR=/private/groups/meyerlab/esteiger

cd ${BASE_DIR}/gvtp_popgen/inversion
python3 cenh3_dnds.py 16
