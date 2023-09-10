# Demultiplexing and doublet detection

This repository contains a pipeline for demultiplexing and doublet detection from Vedolizumab-Predict scRNA-seq data.

'Souporcell_analysis' provides an example of how the functions as a pipeline for assignments of cells to the donor. 

'Hashtag_calling.R' provides the functions for preprocessing scRNA-seq data. The script processes data from CellRanger and integrates Souporcell toolkit results for genotype deconvolution.

## Souporcell Analysis Pipeline 

Clustering mixed-genotype scRNAseq experiments by individual using [souporcell](https://github.com/wheaton5/souporcell)

### Prerequisites
Environment: The pipeline is designed to run on a Slurm cluster.
Software:
souporcell toolkit
PythonPlus
minimap2
SAMtools
freebayes
vartrix
Reference Data: The pipeline uses the GRCh38 human genome reference.

### Pipeline Overview
1. Setup: The script starts by defining Slurm job requirements.
2. Renaming: Using renamer_v2.py, the BAM files are renamed based on the barcodes provided.3
3. Remapping: The renamed sequences are remapped to the reference genome using minimap2.4
4. Retagging: The SAM file from the remapping step is retagged using retag_v2.py.5
5. Sorting & Indexing: The retagged BAM file is sorted and indexed using SAMtools.6
6. Variant Calling: Variants are called using freebayes.7.
7. Allele Counting: Allele counts for the variants are generated using vartrix.8.
8. Clustering: Cells are clustered based on their genotypes using Souporcell.9.
9. Doublet Detection: Doublets (or cells resulting from two cells being sequenced as one) are detected using trouble.

## Hashtag calling Pipeline 

The vignette of demultiplexing with hashtag oligos (HTOs) in [Seurat](https://satijalab.org/seurat/articles/hashing_vignette)


### Prerequisites
R Environment: Ensure you have R installed along with the necessary packages.

Data:
The output from CellRanger (UMI count matrix).
Results from Souporcell for genotype deconvolution.

### Script Overview
1. Read CellRanger Output: The script starts by reading the UMI count matrix using Read10X function.
2. Data Initialization: A Seurat object is created and subsetting is done to keep cells with more than 200 RNA counts.
3. HTO Demultiplexing: Cells are demultiplexed based on HTO data using the HTODemux function.
4. Integrate Souporcell Results: The results from the Souporcell toolkit are integrated to determine the genotype of each cell.
5. QC and Preprocessing:
6. Mitochondrial genes are calculated and visualized.
7. Data is filtered based on QC metrics.
8. Data normalization, scaling, and variable feature selection is performed.
9. Dimensionality Reduction and Clustering:
10. PCA is computed.
11. UMAP is used for dimensionality reduction.
12. Clusters are identified using shared nearest neighbor (SNN) modularity optimization.
13. Visualization: Visualization of clusters using DimPlot and gene expression visualization using FeaturePlot.
14. Save Processed Data: The processed Seurat object is saved for further analysis.
