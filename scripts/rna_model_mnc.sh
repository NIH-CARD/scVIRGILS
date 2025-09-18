#!/bin/bash

# This is so kludgy, but this the input files just need to be based through as arguments
input_file=$1
sample_key=$2
merged_rna_anndata=$3

# Load module
module load singularity/4.1.5
# Run 
singularity run --nv --bind "$PWD" envs/single_cell_gpu.sif python scripts/rna_model_mnc.py "${input_file}" "${sample_key}" "${merged_rna_anndata}"
