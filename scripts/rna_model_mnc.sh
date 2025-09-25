#!/bin/bash
set -euo pipefail

# This is so kludgy, but this the input files just need to be based through as arguments
input_file=$1
seq_batch_key=$2
merged_rna_anndata=$3


# Run 
python scripts/rna_model_mnc.py "${input_file}" "${seq_batch_key}" "${merged_rna_anndata}"
