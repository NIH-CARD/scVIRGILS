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

cell_cycle_genes = work_dir+'/input/lab_cell_cycle_genes.txt'

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
        merged_rna_anndata_unfiltered = work_dir+'/output/01_merged_anndata_rna.h5ad',
        rna_anndata=expand(
            work_dir+'/output/01_{sample}_anndata_object_rna.h5ad', 
            sample=samples
        )

rule preprocess:
    input:
        metadata_table=metadata_table,
        cell_cycle_genes=cell_cycle_genes,
        rna_anndata=lambda wildcards: data_dir+f'data/li_2023/CELLRANGER/{sample_map[wildcards.sample]}/filtered_feature_bc_matrix.h5'
    output:
        rna_anndata = work_dir+'/output/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['single_cell_transcriptomics']
    params:
        sample='{sample}'
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'/scripts/preprocess.py'

rule merge_unfiltered:
    input:
        rna_anndata=expand(
            work_dir+'/output/01_{sample}_anndata_object_rna.h5ad', 
            zip,
            sample=samples
            )
    output:
        merged_rna_anndata_unfiltered = work_dir+'/output/01_merged_anndata_rna.h5ad'
    singularity:
        envs['single_cell_transcriptomics']
    params:
        samples=samples
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'

rule plot_qc_rna:
    input:
        merged_rna_anndata_unfiltered = work_dir+'/output/01_merged_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    resources:
        runtime=960, mem_mb=500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/plot_qc_metrics.py'
