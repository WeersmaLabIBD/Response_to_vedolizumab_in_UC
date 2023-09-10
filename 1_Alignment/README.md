# CellRanger analysis

The 'Cellranger_analysis' contains codes for processing single-cell RNA sequencing (scRNA-seq) data using cellranger version 3.1.0.

## Prerequisites
Environment: The pipeline is designed to run on a Slurm cluster.
Software: cellranger version 3.1.0.
Reference Data: The pipeline uses the GRCh38 human genome reference (version 3.0.0).

## Codes Overview
Setup: The script starts by defining Slurm job requirements such as job name, output and error paths, runtime, memory, CPUs, and temporary storage.
Temporary Directory Creation: The script creates a temporary directory for cellranger operations.
Cell Ranger Analysis: The script runs the cellranger count command using the specified parameters. This command performs demultiplexing, alignment, filtering, and counting.
Moving Results: After cellranger completes its operations, the results are moved from the temporary directory to a final output directory.

## Usage
1. Clone this repository to your local machine.
2. Modify the Slurm directives (#SBATCH lines) at the beginning of the script to fit your cluster's requirements.
3. Ensure the paths for output, error, transcriptome, libraries, and feature references are correctly set in the script.4
4. Submit the script to your Slurm cluster using the sbatch command.
