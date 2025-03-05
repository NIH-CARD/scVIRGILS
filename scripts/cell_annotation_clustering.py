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
sc_marker_genes = snakemake.input.sc_marker_genes

# Step 1: Read in marker genes
print("Reading in marker genes...")

marker_genes_df = pd.read_csv(sc_marker_genes)

print(marker_genes_df)

# Step 3: Read in filtered adata frame
print("Reading in adata object...")

adata = sc.read_h5ad(sc_anndata)

print(adata)