import os
import pandas as pd

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

# Define the data directory, explicitly -- Using Li et al. frontal cortex as a dummy dataset
data_dir = '/data/CARD_MPU/' 

# Define the working directory, explicitly as the directory of this pipeline
work_dir = os.getcwd()

# Number of threads to use when running the rules
num_workers = 8

# Define where the metadata exists for each sample
metadata_table = work_dir+'/input/example_metadata.csv'

# Define where celltypes/cell marker gene information exists
gene_markers_file = work_dir+'/input/example_marker_genes.csv'

# Keys for samples (used in aggregation)
preprocess_key = 'preprocess_id'
sample_key = 'participant_id'

# Read in metadata
df = pd.read_csv(metadata_table)
print("Metadata Columns:", df.columns)  # Debug: Print available column names

# Validate keys before using them
if preprocess_key not in df.columns or sample_key not in df.columns:
    raise KeyError(f"Missing expected columns in metadata: {preprocess_key}, {sample_key}")

preprocess_id = df[preprocess_key].tolist()
samples = df[sample_key].tolist()
sample_map = dict(zip(samples, preprocess_id))  # Map sample → preprocess_id

"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers (downloaded from Quay.io in snakemake.sh)
envs = {
    'single_cell_transcriptomics': 'envs/single_cell_gpu.sif'
}

# Uncomment and correct rule all
rule all:
    input:
        rna_anndata=expand(
            work_dir+'output/01_{sample}_anndata_object_rna.h5ad', 
            sample=samples
        )

rule preprocess:
    input:
        metadata_table=metadata_table,
        rna_anndata=lambda wildcards: data_dir+f'data/li_2023/CELLRANGER/{sample_map[wildcards.sample]}/filtered_feature_bc_matrix.h5'
    output:
        rna_anndata = work_dir+'output/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['single_cell_transcriptomics']
    params:
        sample='{sample}'
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'/scripts/preprocess.py'
