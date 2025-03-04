# Single-cell RNA transcriptomic pipeline

This workflow serves as a generalized pipeline for analyzing single-cell transcriptomic data. This pipeline has been created using a frontal cortex dataset (___enter data reference here____) and is intended to explore gene expression differences between FTD, AD, and control samples across different cell types.  

fastqs --> Align -- BAM (cellranger) --> Remove doublets --> unfiltered h5ad --> filtered h5ad --> cell clustering/annotation --> cell extraction --> pseudobulking (Seurat) --> Normalization (Seurat) --> Batch correction (Seurat, CombatSeq, CCA, etc) --> DEGs (DEseq2), GO (gprofiler2), cell-type proportions, etc.

## Required Input

Metadata file in .csv format, example in input/example_metadata.csv. A minimal metadata file should include:

- 






- Fastq files
- Marker genes
