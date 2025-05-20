import scanpy as sc

# Read in data
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Split out only the celltype AnnData, and reduce the size of the object
adata = adata[adata.obs['cell_type'] == snakemake.params.cell_type].copy()

# Write the rna data out
adata.write_h5ad(snakemake.output.celltype_rna, compression='gzip')