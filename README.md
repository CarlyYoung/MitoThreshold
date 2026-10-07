# MitoThreshold
Description: An R pipeline for evaluating mitochondrial RNA-content thresholds in single-cell RNA-seq data

The pipeline processes unfiltered, unintegrated Seurat data, performs split-aware SCTransform preprocessing and DoubletFinder detection, and evaluates response trajectories across a range of mitochondrial-content thresholds. It reports data-driven transition points for multiple technical and biological metrics.

## Key analyses
- Cell retention and library complexity
- Mean gene-expression variance
- PCA variance
- Housekeeping, apoptosis, mitochondrial and lineage-marker expression
- Cell-level and pseudobulk AUCell pathway activity
- Automated detection of mitochondrial-threshold transition points

## Input
The pipeline accepts:
- A 10x feature-barcode `.h5` file
- An `.rds` file containing a Seurat object
- An `.RData` or `.rda` file containing a Seurat object

Raw RNA counts are required. Sample, split and optional biological-group metadata are configured in the script.

## Requirements
The pipeline uses R and the following packages:
`Seurat`, `DoubletFinder`, `dplyr`, `purrr`, `tidyr`, `stringr`, `tibble`, `ggplot2`, `Matrix`, `slider`, `readr`, `future`, `AUCell` and `msigdbr`.

The `hdf5r` package is also required when importing 10x `.h5` files.

## Usage
Edit the `USER CONFIGURATION` section in `MitoThreshold_Reproducible_Pipeline.R`, then run:

```bash
Rscript MitoThreshold_Reproducible_Pipeline.R
```

Alternatively, run the script from a fresh RStudio session:

```r
source("MitoThreshold_Reproducible_Pipeline.R")
```

Results are written to separate preprocessing and mitochondrial-threshold analysis directories. Existing non-empty output directories are never overwritten.

## Citation
Citation details for the associated manuscript will be added upon publication.
