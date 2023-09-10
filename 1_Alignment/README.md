# CellRanger analysis

The starting point for the analysis is single-cell RNA-seq (scRNA-seq) data generated with the Chromium 10X and mapped using [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).

The 'Cellranger_analysis' contains codes for processing scRNA-seq data using cellranger version 3.1.0.

## Prerequisites
Environment: The pipeline is designed to run on a Slurm cluster.
Software: cellranger version 3.1.0.
Reference Data: The pipeline uses the GRCh38 human genome reference (version 3.1.0).

## Pipeline Overview
1. Setup: The script starts by defining Slurm job requirements such as job name, output and error paths, runtime, memory, CPUs, and temporary storage.
2. Temporary Directory Creation: The script creates a temporary directory for cellranger operations.
3. Cell Ranger Analysis: The script runs the cellranger count command using the specified parameters. This command performs demultiplexing, alignment, filtering, and counting.
4. Moving Results: After cellranger completes its operations, the results are moved from the temporary directory to a final output directory.
