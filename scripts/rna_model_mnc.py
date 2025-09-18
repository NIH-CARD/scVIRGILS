import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np
import anndata as ad
import snapatac2 as snap
import scanpy.external as sce

# Read in AnnData atlas object
adata = ad.read_h5ad(sys.argv[1])

# Plot UMAPs
# Compute overall UMAP
sc.pp.neighbors(adata, n_neighbors=10)  # Compute neighbors
sc.tl.umap(adata)  # Run UMAP

# Check available representations
snap.pp.select_features(adata, n_features=5000, n_jobs=30)

# Spectral MDS analysis - reduces the dimensionality of the dataset by finding the simple representation of the data
snap.tl.spectral(adata) # this adds adata.obsm["X_spectral"]

# Batch correction - This is a modified MNN correct algorithm based on cluster centroid (relative to mnn_correct())
snap.pp.mnc_correct(adata, batch=sys.argv[2]) #This adds adata_cluster5.obsm["X_spectral_mnn"]

# Perform k-nearest neighbors - This is done after batch correction to match similar cells based on corrected accessibility signals. The data is more comparable for this step after MNC corrections
snap.pp.knn(adata, use_rep='X_spectral_mnn') #This runs and replaces adata_cluster5.obsm["X_spectral_mnn"]

# Cluster 
snap.tl.leiden(adata)
# Calculate UMAP
snap.tl.umap(adata, use_rep='X_spectral_mnn', key_added='umap_mnn') # This adds adata.obsm["X_umap_mnn"]

# Plot UMAPs
# Compute overall UMAP
sc.pp.neighbors(adata, n_neighbors=10)  # Compute neighbors
sc.tl.umap(adata)  # Run UMAP

# Save the anndata object
adata.write_h5ad(sys.argv[4], compression='gzip')