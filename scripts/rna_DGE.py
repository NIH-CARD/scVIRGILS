import anndata as ad
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import decoupler as dc
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

seq_batch_key = snakemake.params.seq_batch_key

# Open rna
adata = sc.read_h5ad(snakemake.input.rna_anndata)

# Convert SampleID to categorical and then to codes (integers)
adata.obs[seq_batch_key] = pd.Categorical(adata.obs[seq_batch_key]).codes

adata = adata[adata.obs['cell_type'] == snakemake.params.cell_type].copy()

# Read in parameters
disease_name = snakemake.params.disease
control_name = snakemake.params.control
disease_param = snakemake.params.disease_param

# Get pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col=snakemake.params.sample_key,
    groups_col=disease_param,
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=10
)

# Store raw counts in layers
pdata.layers['counts'] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)

# Return raw counts to X
dc.swap_layer(pdata, 'counts', X_layer_key=None, inplace=True)

# Abbreviate diagnosis to avoid space syntax error
pdata.obs['comparison'] = pdata.obs[disease_param]

dc.get_metadata_associations(
    pdata,
    obs_keys = ['comparison', 'psbulk_n_cells', 'psbulk_counts'],  # Metadata columns to associate to PCs
    obsm_key='X_pca',  # Where the PCs are stored
    uns_key='pca_anova',  # Where the results are stored
    inplace=True,
)

# Export pseudobulk
# CSV pseudobulk
adata_df = pd.DataFrame(pdata.X)
sample_cell = pdata.obs[[snakemake.params.sample_key, 'cell_type', disease_param]]
adata_df.columns = pdata.var_names.to_list()
adata_df.index = sample_cell.index
adata_df = pd.merge(left=sample_cell, right=adata_df, left_index=True, right_index=True)
adata_df.to_csv(snakemake.output.celltype_pseudobulk, index=False)

#pdata.write_h5ad(snakemake.output.celltype_pseudobulk)

pdata_genes = dc.filter_by_expr(
    pdata, 
    group='comparison', 
    min_count=10, 
    min_total_count=15
    )

# Subset valuable genes
pdata = pdata[:, pdata_genes].copy()

# Determine the number of cpus to use
inference = DefaultInference(n_cpus=1)

dds = DeseqDataSet(
    adata=pdata,
    design_factors=['comparison', snakemake.params.seq_batch_key],
    refit_cooks=True,
    inference=inference,
)

# Compute log-fold changes
dds.deseq2()

# Extract contrast between control and disease states
stat_res = DeseqStats(
    dds,
    contrast=['comparison', disease_name, control_name],
    inference=inference,
)

# Compute Wald test
stat_res.summary()

# Extract results
DGE_results_df = stat_res.results_df
DGE_results_df['-log10_padj'] = -np.log10(DGE_results_df['padj'])
DGE_results_df.to_csv(snakemake.output.output_DGE_data)

# Plot 
dc.plot_volcano_df(
    DGE_results_df,
    x='log2FoldChange',
    y='padj',
    top=20,
    lFCs_thr=1,
    sign_thr=1e-2,
    figsize=(4, 4),
    dpi=600,
    save=snakemake.output.output_figure
)
