#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 3) {
  cat("Usage:\n",
      "  Rscript h5ad_to_seurat_seuratdisk.R <input.h5ad> [output.rds] [assay=RNA]\n",
      sep = "")
  quit(status = 1)
}

in_h5ad <- normalizePath(args[1])
if (!file.exists(in_h5ad)) stop("Input not found: ", in_h5ad)
out_rds <- if (length(args) >= 2) args[2] else sub("\\.h5ad$", ".rds", basename(in_h5ad), ignore.case = TRUE)
assay   <- if (length(args) >= 3) args[3] else "RNA"

# Work in the input's directory so Convert() writes the .h5seurat next to it
wd_old <- getwd()
on.exit(setwd(wd_old), add = TRUE)
setwd(dirname(in_h5ad))
h5s_out <- sub("\\.h5ad$", ".h5seurat", basename(in_h5ad), ignore.case = TRUE)

message("Converting to h5seurat …")
Convert(basename(in_h5ad), dest = "h5seurat", assay = assay, overwrite = TRUE, verbose = TRUE)

