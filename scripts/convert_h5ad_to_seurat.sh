#!/bin/bash
set -euo pipefail

# Resolve repository root as the parent of this script (adjust if your layout differs)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ATLAS_DIR="$REPO_ROOT/atlas"

# Ensure required files exist up-front (clear error if not)
if [[ ! -f "$ATLAS_DIR/03_modeled_anndata_rna.h5ad" ]]; then
  echo "ERROR: Missing $ATLAS_DIR/03_modeled_anndata_rna.h5ad" >&2
  exit 1
fi
if [[ ! -f "$REPO_ROOT/scripts/convert_h5ad_to_seurat.R" ]]; then
  echo "ERROR: Missing $REPO_ROOT/scripts/convert_h5ad_to_seurat.R" >&2
  exit 1
fi
if [[ ! -f "$REPO_ROOT/envs/r_container.sif" ]]; then
  echo "ERROR: Missing $REPO_ROOT/envs/r_container.sif" >&2
  exit 1
fi

# Load modules (works when run via bash -l/-lc or in environments with modulecmd)
if command -v module >/dev/null 2>&1; then
  module load singularity
fi

# Bind the *host* atlas directory into the container at /atlas
# Then call the R script with explicit paths so CWD doesn’t matter.
singularity exec \
  --bind "$REPO_ROOT:$REPO_ROOT" \
  "$REPO_ROOT/envs/r_container.sif" \
  Rscript "$REPO_ROOT/scripts/convert_h5ad_to_seurat.R" "$REPO_ROOT/atlas/03_modeled_anndata_rna.h5ad"
