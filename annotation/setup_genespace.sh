#!/bin/bash
# One-time setup. Builds conda env genespace: r-base=4.3, MCScanX, GENESPACE's
# Bioconductor/CRAN dependencies, then GENESPACE 1.3.1 from GitHub.
# GENESPACE is fed a precomputed OrthoFinder run (run_orthofinder_genespace.sh)
# rather than calling it.
# Safe to re-run.
#
# Usage: bash setup_genespace.sh

set -euo pipefail

CONDA_BASE=/private/home/esteiger/miniconda3
ENV_NAME=genespace

source "${CONDA_BASE}/etc/profile.d/conda.sh"

if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo "Conda env '${ENV_NAME}' already exists — skipping creation."
else
    echo "Creating conda env '${ENV_NAME}'..."
    conda create -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
        r-base=4.3 r-data.table r-dbscan bioconductor-biostrings r-igraph \
        r-r.utils r-ggplot2 r-devtools r-ape bioconductor-rtracklayer \
        bioconductor-genomicranges bioconductor-iranges \
        bioconductor-biocgenerics bioconductor-genomeinfodb r-gridextra \
        mcscanx
fi

conda activate "${ENV_NAME}"

# path2mcscanx must be the directory containing MCScanX_h, i.e. the env's bin/
if [ ! -x "${CONDA_PREFIX}/bin/MCScanX_h" ]; then
    echo "ERROR: MCScanX_h not found at ${CONDA_PREFIX}/bin/MCScanX_h after install" >&2
    exit 1
fi
echo "MCScanX_h found: ${CONDA_PREFIX}/bin/MCScanX_h"

"${CONDA_PREFIX}/bin/Rscript" -e '
if (!requireNamespace("GENESPACE", quietly = TRUE)) {
  devtools::install_github("jtlovell/GENESPACE", upgrade = "never")
} else {
  cat("GENESPACE already installed (version", as.character(packageVersion("GENESPACE")), ") - skipping.\n")
}
'

echo "Setup complete. Env: ${ENV_NAME}"
echo "path2mcscanx: ${CONDA_PREFIX}/bin"
echo "Next: run_orthofinder_genespace.sh must have completed, then run setup_genespace_run.R"
