import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np

# Keep consistent font sizes

SMALL_SIZE = 6
MEDIUM_SIZE = 8
BIGGER_SIZE = 10

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE) 

# Convert to metric
cm = 1/2.54

# Set figure context
sns.set_context('paper')

# Save the AnnData object
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata_unfiltered) 

for sample in adata.obs['sample'].drop_duplicates().to_list():

    # Make plot directory
    try:
        os.mkdir(f'plots/{sample}')
    except FileExistsError:
        print('Already there')


    """Plot percent mitochondria 
    """
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
    fig.suptitle(f' Sample {sample} ', fontsize=BIGGER_SIZE)

    # Violin plot in the first panel
    sc.pl.violin(adata[adata.obs['sample'] == sample], ['pct_counts_mt'], jitter=0.5, ax=ax[0], show=False)
    ax[0].set_ylabel('percent')
    ax[0].set_xticks('')
    ax[0].set_xlim(-.75, .75)
    ax[0].plot([-.5, .5], [20, 20], '--r')
    ax[0].set_title('Percent mitochondria per cell')

    # Histogram of values in the second panel
    y, x, _ = ax[1].hist(
        adata[adata.obs['sample'] == sample].obs['pct_counts_mt'], 
        bins=int(np.sqrt(adata[adata.obs['sample'] == sample].n_obs))
        )
    ax[1].set_yscale("log")
    ax[1].set_xlabel('percent')
    ax[1].set_ylabel('number of cells')
    ax[1].plot([15, 15], [1, y.max()], '--r')
    ax[1].set_ylim(0, y.max())
    ax[1].set_title('Percent mitochondria per cell')
    plt.savefig(f'plots/{sample}/mito_pct.png', dpi=300)


    """Plot percent ribosome"""
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(f' Sample {sample} ', fontsize=BIGGER_SIZE)

    # Violin plot in the first panel
    sc.pl.violin(adata[adata.obs['sample'] == sample], ['pct_counts_rb'], jitter=0.5, ax=ax[0], show=False)
    ax[0].set_ylabel('percent')
    ax[0].set_xlim(-.75, .75)
    ax[0].plot([-.5, .5], [15, 15], '--r')
    ax[0].set_title('Percent ribosome genes per cell')

    # Histogram of values in the second panel
    y, x, _ = ax[1].hist(
        adata[adata.obs['sample'] == sample].obs['pct_counts_rb'], 
        bins=int(np.sqrt(adata[adata.obs['sample'] == sample].n_obs))
        )
    ax[1].plot([10, 10], [1, y.max()], '--r')
    ax[1].set_ylim(0, y.max())
    ax[1].set_yscale("log")
    ax[1].set_xlabel('percent')
    ax[1].set_ylabel('number of cells')
    ax[1].set_title('Percent ribosome genes per cell')
    plt.savefig(f'plots/{sample}/ribo_pct.png', dpi=300)


    """Plot number of genes per cell"""
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(f' Sample {sample} ', fontsize=BIGGER_SIZE)

    # Violin plot in the first panel
    sc.pl.violin(adata[adata.obs['sample'] == sample], ['n_genes_by_counts'], jitter=0.5, ax=ax[0], show=False)
    ax[0].plot([-.5, .5], [500, 500], '--r')
    ax[0].set_ylabel('total counts')
    ax[0].set_xlim(-.75, .75)
    ax[0].plot([-.5, .5], [500, 500], '--r')
    ax[0].set_title('Number of genes per cell')

    # Histogram of values in the second panel
    y, x, _ = ax[1].hist(
        adata[adata.obs['sample'] == sample].obs['n_genes_by_counts'], 
        bins=int(np.sqrt(adata[adata.obs['sample'] == sample].n_obs))
        )
    ax[1].plot([500, 500], [1, y.max()], '--r')
    ax[1].set_ylim(0, y.max())
    ax[1].set_xlabel('total counts')
    ax[1].set_ylabel('number of cells')
    ax[1].set_title('Number of genes per cell')
    plt.savefig(f'plots/{sample}/num_genes_per_cell.png', dpi=300)


    """Plot the scrublet values"""
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(f' Sample {sample} ', fontsize=BIGGER_SIZE)

    # Violin plot in the first panel
    sc.pl.violin(adata[adata.obs['sample'] == sample], ['doublet_score'], jitter=0.5, ax=ax[0], show=False)
    ax[0].plot([-.5, .5], [0.25, 0.25], '--r')
    ax[0].set_ylabel('droplet score')
    ax[0].set_xlim(-.75, .75)
    ax[0].plot([-.5, .5], [0.25, 0.25], '--r')
    ax[0].set_title('Doublet score per cell')

    
    y, x, _ = ax[1].hist(
        adata[adata.obs['sample'] == sample].obs['doublet_score'], 
        bins=int(np.sqrt(adata[adata.obs['sample'] == sample].n_obs))
        )
    ax[1].plot([0.25, 0.25], [1, y.max()], '--r')
    ax[1].plot([0.25, 0.25], [1, y.max()], '--r')
    ax[1].set_ylim(0, y.max())
    ax[1].set_xlabel('droplet score')
    ax[1].set_ylabel('number of droplets')
    ax[1].set_title('Doublet score per cell')
    plt.savefig(f'plots/{sample}/scrublet_score_per_cell.png', dpi=300)

    sc.pl.scatter(
        adata[adata.obs['sample'] == sample], 
        "total_counts", 
        "n_genes_by_counts", 
        color="pct_counts_mt",
        show=False
        )
    plt.savefig(f'plots/{sample}/num_gene_counts_total.png', dpi=300)

    # Plot UMAP of values
    sc.pl.umap(
        adata[adata.obs['sample'] == sample],
        color=['pct_counts_mt', 'pct_counts_rb', 'doublet_score', 'total_counts', 'n_genes_by_counts'],
        size=2,
        ncols=3,
        cmap='viridis',
        show=False
    )
    plt.savefig(f'plots/{sample}/qc_umap.png', dpi=300)