#!/bin/bash
# One-time setup. Builds conda env orthofinder:
#   conda create -n orthofinder -c bioconda -c conda-forge orthofinder=3.1.5
# Safe to re-run; run on a login node.
#
# Usage: bash setup_orthofinder.sh

set -euo pipefail

CONDA_BASE=/private/home/esteiger/miniconda3
ENV_NAME=orthofinder

source "${CONDA_BASE}/etc/profile.d/conda.sh"

if conda env list | grep -qE "^${ENV_NAME}\s"; then
    echo "Conda env '${ENV_NAME}' already exists — skipping creation."
else
    echo "Creating conda env '${ENV_NAME}'..."
    conda create -y -n "${ENV_NAME}" -c bioconda -c conda-forge orthofinder
fi

# the aster package's activate.d script fails under set -u without this
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

conda activate "${ENV_NAME}"

# call the env's python explicitly; /bin precedes the env on PATH here
"${CONDA_PREFIX}/bin/python3" "${CONDA_PREFIX}/bin/orthofinder" -h | head -5

# OrthoFinder 3.1.5's --assign path calls np.chararray, removed in numpy 2;
# patch that line rather than downgrade numpy
SAMPLE_GENES="${CONDA_PREFIX}/bin/src/orthofinder/orthogroups/sample_genes.py"
if grep -q "np.chararray((m, len(I)))" "${SAMPLE_GENES}"; then
    echo "Patching np.chararray -> np.empty(..., dtype='<U1') in ${SAMPLE_GENES} (NumPy 2 compat)..."
    cp -n "${SAMPLE_GENES}" "${SAMPLE_GENES}.orig-pre-numpy2-patch"
    sed -i "s/z = np.chararray((m, len(I)))/z = np.empty((m, len(I)), dtype='<U1')/" "${SAMPLE_GENES}"
else
    echo "sample_genes.py already patched (or no longer uses np.chararray) — skipping."
fi
echo "Setup complete. Env: ${ENV_NAME}"
echo "Next: for any species not yet processed, run"
echo "  prep_ensembl_plant_proteome.sh (Ensembl Plants) or"
echo "  prep_ncbi_datasets_proteome.sh (NCBI Datasets, e.g. Helianthus anomalus)"
echo "  then symlink its .pep.longest.fa into:"
echo "    orthofinder_input/          (full accumulated set — only used for a from-scratch -f rerun)"
echo "    orthofinder_input_assign/   (ONLY the new species not yet folded into a --core run — see"
echo "                                 run_orthofinder.sh for why OrthoFinder 3 uses --core/--assign"
echo "                                 instead of the old -b/-f 'add species' workflow)"
echo "  before: sbatch run_orthofinder.sh"
