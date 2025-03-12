# Single-cell RNA transcriptomic pipeline

This workflow serves as a generalized pipeline for analyzing single-cell transcriptomic data. This pipeline has been created using a frontal cortex dataset (___enter data reference here____) and is intended to explore gene expression differences between FTD, AD, and control samples across different cell types.  

fastqs --> Align -- BAM (cellranger) --> Remove doublets --> unfiltered h5ad --> filtered h5ad --> cell clustering/annotation --> cell extraction --> pseudobulking (Seurat) --> Normalization (Seurat) --> Batch correction (Seurat, CombatSeq, CCA, etc) --> DEGs (DEseq2), GO (gprofiler2), cell-type proportions, etc.

## To get started

Copy this repository to where you will be working with your dtaa. This folder will be where output data is stored, while intermediary files will be stored in a separate folder to be defined by the user.

### Required Input

Metadata file in .csv format, example in input/example_metadata.csv. A minimal metadata file should include:

- 
- Fastq files

Marker genes file in .csv format, example in input/example_marker_genes.csv. A minimal marker gene file sould include:

- official gene symbol
- cell type


Once set up, this complete pipeline can be run with `bash snakemake.sh` in terminal. **Note: This can only be run with an interactive or slurm job.**