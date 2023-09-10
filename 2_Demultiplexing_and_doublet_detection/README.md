# Demultiplexing and doublet detection

This repository contains a pipeline for demultiplexing and doublet detection from Vedolizumab-Predict scRNA-seq data.

'Souporcell_analysis' provides an example of how the functions as a pipeline for assignments of cells to the donor. 

'Hashtag_calling.R' provides the functions for demultiplexing and doublet detection.

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

### Usage
1. Clone this repository to your local machine.
2. Modify the Slurm directives (#SBATCH lines) at the beginning of the script to fit your cluster's requirements.3.
3. Ensure the paths for output, error, transcriptome, libraries, and other references are correctly set in the script.
4. Submit the script to your Slurm cluster using the sbatch command.

