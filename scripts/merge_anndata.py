import anndata as ad
import pandas as pd
import scanpy as sc

# Create dictionary of sample name to file, then filter with relevant diagnosis
sample_loc = dict(zip(snakemake.params.samples, snakemake.input.rna_anndata))
adatas = {key: sc.read_h5ad(sample_loc[key]) for key in snakemake.params.samples}

# Concatenate with names of sample
adata = ad.concat(
    merge='same', index_unique='_', join='outer',
    adatas=adatas
    )

# Write out the unfiltered dataset
adata.write_h5ad(
    filename=snakemake.output.merged_rna_anndata_unfiltered, 
    compression='gzip'
    )