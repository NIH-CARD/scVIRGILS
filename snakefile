import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

# Define the data directory, explicitly -- Using Li et al. frontal cortex as a dummy dataset
data_dir = '/data/CARD_MPU/ftd-ad-atlas/datasets/li_2023/objects/' # may have to 
# Define the working directory, explictly as the directory of this pipeline
work_dir = os.getcwd()

# Number of threads to use when running the rules
num_workers = 8

# Define where the metadata data exists for each sample to be processed
metadata_table = work_dir+'/input/example_metadata.csv'
# Define where celltypes/cell marker gene 
gene_markers_file = work_dir+'/input/example_marker_genes.csv'

"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'single_cell_transcriptomics': 'envs/single_cell_cpu.sif'
    }


rule annotate:
    input:
        merged_rna_anndata = data_dir+'li_filtered.h5ad',
        gene_markers = gene_markers_file
    output:
        merged_rna_anndata = work_dir+'/output/li_annotated.h5ad',
        cell_annotate = work_dir+'/output/rna_cell_annot.csv'
    singularity:
        envs['single_cell_transcriptomics']
    resources:
        runtime=240, mem_mb=500000, slurm_partition='largemem'
    script:
        'scripts/annotate.py'