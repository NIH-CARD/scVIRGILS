import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

# Define the data directory, explicitly -- Using Li et al. frontal cortex as a dummy dataset
data_dir = '/data/CARD_MPU/ftd-ad-atlas/datasets/li_2023/objects/'
# Define the working directory, explictly as the directory of this pipeline
work_dir = os.getcwd()

# Number of threads to use when running the rules
num_workers = 8


"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'single_cell_transcriptomics': 'envs/single_cell_cpu.sif'
    }


rule annotate:
    input:
        sc_anndata = data_dir+'li_filtered.h5ad'
    output:
        sc_annotated_h5ad = data_dir+'li_annotated.h5ad',
        sc_annotated_plot = data_dir+'umap_cell_type.png'
    singularity:
        envs['single_cell_transcriptomics']
    script:
        work_dir+'/scripts/cell_annotation_clustering.py'