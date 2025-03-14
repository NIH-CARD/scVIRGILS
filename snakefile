import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

# Define the data directory, explicitly -- Using Li et al. frontal cortex as a dummy dataset
data_dir = '/data/CARD_MPU/' 
# Define the working directory, explictly as the directory of this pipeline
work_dir = os.getcwd()

# Number of threads to use when running the rules
num_workers = 8

# Define where the metadata data exists for each sample to be processed
metadata_table = work_dir+'/input/example_metadata.csv'
# Define where celltypes/cell marker gene 
gene_markers_file = work_dir+'/input/example_marker_genes.csv'

# Key for samples, required in aggregating while preserving sample info
preprocess_key = 'preprocess_id'
sample_key = 'participant_id'

# Read in the list of batches and samples
#batches = pd.read_csv(metadata_table)['Use_batch'].tolist()
preprocess_id = pd.read_csv(metadata_table)[preprocess_key].tolist()
samples = pd.read_csv(metadata_table)[sample_key].tolist()


"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'single_cell_transcriptomics': 'envs/single_cell_cpu.sif'
    }

rule all:
    input:
        rna_anndata=expand(
            work_dir+'output/01_{samples}_anndata_object_rna.h5ad', 
            zip,
            preprocess_id=preprocess_id,
            sample=samples
            )


rule preprocess:
    input:
        metadata_table=metadata_table,
        rna_anndata = data_dir+'data/li_2023/CELLRANGER/{preprocess_id}/filtered_feature_bc_matrix.h5'
    output:
        rna_anndata = work_dir+'output/01_{samples}_anndata_object_rna.h5ad'
    singularity:
        envs['single_cell_transcriptomics']
    params:
        sample='{sample}'
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'/scripts/preprocess.py'

#rule annotate:
#    input:
#        merged_rna_anndata = data_dir+'li_filtered.h5ad',
#        gene_markers = gene_markers_file
#    output:
#        merged_rna_anndata = work_dir+'/output/li_annotated.h5ad',
#        cell_annotate = work_dir+'/output/rna_cell_annot.csv'
#    singularity:
#        envs['single_cell_transcriptomics']
#    resources:
#        runtime=240, mem_mb=500000, slurm_partition='largemem'
#    script:
#        'scripts/annotate.py'