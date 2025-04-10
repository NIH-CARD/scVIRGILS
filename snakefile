import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

"""File locations"""
data_dir = 'YOUR DIRECTORY HERE' # Define the data directory, explicitly
work_dir = os.getcwd() # Define the working directory, explictly as the directory of this pipeline
metadata_table = work_dir+'/input/example_metadata.csv' # Define where the metadata data exists for each sample to be processed
gene_markers_file = work_dir+'/input/example_marker_genes.csv' # Define where celltypes/cell marker gene 

"""Metadata parameters"""
seq_batch_key = 'Use_batch' # Key for sequencing batch, used for directory search
sample_key = 'Sample' # Key for samples, required in aggregating while preserving sample info
batches = pd.read_csv(metadata_table)[seq_batch_key].tolist() # Read in the list of batches and samples
samples = pd.read_csv(metadata_table)[sample_key].tolist()
disease_param = 'Primary Diagnosis' # Name of the disease parameter
control = 'control' # Define disease states
diseases = ['PD', 'DLB'] # Disease states to compare, keep as list of strings, unnecessary 
cell_types = pd.read_csv(gene_markers_file)['cell type'] # Define the cell types to look for, from gene marker file

"""Quality control thresholds"""
mito_percent_thresh = 15 # Maximum percent of genes in a cell that can be mitochondrial
ribo_percent_thresh = 10 # Maximum percent of genes in a cell that can be ribosomal
doublet_thresh = 0.15 # Maximum doublet score for a cell, computed by scrublet
min_genes_per_cell = 250 # Minimum number of unique genes in a cell

"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Singularity containers to be downloaded from Quay.io, done in snakemake.sh
envs = {
    'snapatac2': 'envs/snapatac2.sif',
    'singlecell': 'envs/single_cell_gpu.sif',
    'scenicplus': 'envs/scenicplus.sif',
    'decoupler': 'envs/decoupler.sif'
    }

rule all:
    input:
        genes_by_counts = work_dir+'figures/QC_genes_by_counts.png
# Uncomment to view QC data
"""genes_by_counts = work_dir+'figures/QC_genes_by_counts.png'"""
# Uncomment when you have verified QC metrics
"""rna_anndata=expand(
            data_dir+'/{sample}/03_{sample}_anndata_filtered_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            ),"""
# Uncomment when you want to model rna data
"""merged_rna_anndata = work_dir+'/atlas/05_annotated_anndata_rna.h5ad'"""
# Uncomment when you want to run DGE/DAR analysis
"""output_DGE_data = expand(
    work_dir + '/data/significant_genes/rna/rna_{cell_type}_{disease}_DGE.csv',
    cell_type = cell_types,
    disease = diseases
    )
"""
rule cellbender:
    input:
        rna_anndata =data_dir+'/{sample}/raw_feature_bc_matrix.h5',
        cwd = data_dir+'/{sample}/'
    output:
        rna_anndata = data_dir+'/{sample}/cellbender_gex_counts_filtered.h5'
    params:
        sample='{sample}'
    resources:
        runtime=2880, mem_mb=300000, gpu=1, gpu_model='v100x'
    shell:
        work_dir+'/scripts/cellbender_array.sh {input.rna_anndata} {input.cwd} {output.rna_anndata}'

rule rna_preprocess:
    input:
        metadata_table=metadata_table,
        rna_anndata = data_dir+'/{sample}/cellbender_gex_counts_filtered.h5'
    output:
        rna_anndata = data_dir+'/{sample}/01_{sample}_anndata_object_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        sample='{sample}',
        sample_key = sample_key
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'/scripts/rna_preprocess.py'

rule merge_unfiltered:
    input:
        rna_anndata=expand(
            data_dir+'/{sample}/01_{sample}_anndata_object_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_rna_anndata = work_dir+'/atlas/01_merged_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        samples=samples
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'

rule plot_qc_rna:
    input:
        merged_rna_anndata = work_dir+'/atlas/01_merged_anndata_rna.h5ad'
    output:
        mito_figure = work_dir+'/figures/QC_mito_pct.png',
        ribo_figure = work_dir+'/figures/QC_ribo_pct.png',
        gene_counts_figure = work_dir+'/figures/QC_gene_counts.png',
        doublet_figure = work_dir+'/figures/QC_doublet.png',
        genes_by_counts = work_dir+'figures/QC_genes_by_counts.png'
    singularity:
        envs['singlecell']
    resources:
        runtime=960, mem_mb=500000, disk_mb=10000, slurm_partition='largemem' 
    params:
        mito_percent_thresh = mito_percent_thresh,
        doublet_thresh = doublet_thresh,
        min_genes_per_cell = min_genes_per_cell,
        ribo_percent_thresh = ribo_percent_thresh,
        sample_key=sample_key,
        
    script:
        work_dir+'/scripts/plot_qc_metrics.py'

rule filter_rna:
    input:        
        rna_anndata = data_dir+'/{sample}/01_{sample}_anndata_object_rna.h5ad'
    output:
        rna_anndata = data_dir+'/{sample}/02_{sample}_anndata_filtered_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        mito_percent_thresh = mito_percent_thresh,
        doublet_thresh = doublet_thresh,
        min_genes_per_cell = min_genes_per_cell,
        ribo_percent_thresh = ribo_percent_thresh
    resources:
        runtime=120, mem_mb=100000, disk_mb=10000, slurm_partition='quick' 
    script: 
        work_dir+'/scripts/rna_filter.py'

rule merge_filtered_rna:
    input:
        rna_anndata=expand(
            data_dir+'/{sample}/02_{sample}_anndata_filtered_rna.h5ad', 
            zip,
            batch=batches,
            sample=samples
            )
    output:
        merged_rna_anndata = work_dir+'/atlas/02_filtered_anndata_rna.h5ad'
    singularity:
        envs['singlecell']
    params:
        samples=samples
    resources:
        runtime=120, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'/scripts/merge_anndata.py'

rule rna_model:
    input:
        merged_rna_anndata = work_dir+'/atlas/02_filtered_anndata_rna.h5ad'
    output:
        merged_rna_anndata = work_dir+'/atlas/03_modeled_anndata_rna.h5ad',
        model_history = work_dir+'/data/model_elbo/rna_model_history.csv'
    params:
        model = work_dir+'/data/models/rna/',
        sample_key = sample_key
    threads:
        64
    resources:
        runtime=2880, mem_mb=300000, gpu=2, gpu_model='v100x'
    shell:
        'scripts/rna_model.sh {input.merged_rna_anndata} {params.sample_key} {output.model_history} {output.merged_rna_anndata} {params.model}'

rule annotate:
    input:
        merged_rna_anndata = work_dir+'/atlas/03_modeled_anndata_rna.h5ad',
        gene_markers = gene_markers_file
    output:
        merged_rna_anndata = work_dir+'/atlas/04_annotated_anndata_rna.h5ad',
        cell_annotate = work_dir+'/data/rna_cell_annot.csv'
    singularity:
        envs['singlecell']
    params:
        seq_batch_key=seq_batch_key
    resources:
        runtime=240, mem_mb=500000, slurm_partition='largemem'
    script:
        'scripts/annotate.py'

rule export_celltypes:
    input:
        merged_rna_anndata = work_dir+'/atlas/04_annotated_anndata_rna.h5mu'
    output:
        celltype_rna = work_dir+'/data/celltypes/{cell_type}/rna.h5ad'
    params:
        cell_type = lambda wildcards, output: output[0].split('/')[-2]
    singularity:
        envs['singlecell']
    threads:
        8
    resources:
        runtime=120, mem_mb=300000
    script:
        'scripts/export_celltype.py'

rule DGE:
    input:
        rna_anndata = work_dir + '/atlas/05_annotated_anndata_rna.h5ad'
    output:
        output_DGE_data = work_dir + '/data/significant_genes/rna/rna_{cell_type}_{disease}_DGE.csv',
        output_figure = work_dir + 'figures/{cell_type}/rna_{cell_type}_{disease}_DAR.png',
        celltype_pseudobulk = work_dir+'/data/celltypes/{cell_type}/rna_{cell_type}_{disease}_pseudobulk.h5ad'
    params:
        disease_param = disease_param,
        control = control,
        disease = lambda wildcards, output: output[0].split("_")[-2],
        cell_type = lambda wildcards, output: output[0].split("_")[-3]
    singularity:
        envs['decoupler']
    threads:
        64
    resources:
        runtime=1440, disk_mb=200000, mem_mb=200000
    script:
        'scripts/rna_DGE.py'