#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp
import anndata as ad
import numpy as np
import os

sc_anndata = snakemake.input.sc_anndata

# Step 1: Define the marker genes and their associated cell types - Only using those gene sthat have been tested in FCX from spreadsheet
print("defining marker genes...")

marker_genes_dict = {
    "Neuron": ["GRIN1", "GRIN2A", "GAD2", "NEUROD6", "THEMIS"],
    "Astrocyte": ["COL5A3", "GFAP", "SLC1A2", "AQP4", "PDGFRB"],
    "Oligodendrocyte": ["CLDN11", "MOBP", "ST18", "MOG"],
    "Oligodendrocyte precursor": ["LHFPL3", "MEGF11", "PCDH15", "VCAN"],
    "Immune Cell": ["PTPRC", "CSF1R", "CD4", "ITGAM", "P2RY12", "MS4A1", "CD8B", "CD19"],
    "Brain Vascular Cell": ["VCAM1", "CD34", "FLT1", "NOS3", "LCN6", "ANPEP", "FOXC1"],
}

# Step 2: Convert the marker genes dictionary to a DataFrame
print("Converting marker genes dictionary to dataframe...")

marker_genes = [(gene, cell_type) for cell_type, genes in marker_genes_dict.items() for gene in genes]
marker_genes_df = pd.DataFrame(marker_genes, columns=["gene_name", "cell_type"])

print(marker_genes_df)

# Step 3: Read in filtered adata frame
print("Reading in adata object...")

adata = sc.read_h5ad(sc_anndata)

print(adata)