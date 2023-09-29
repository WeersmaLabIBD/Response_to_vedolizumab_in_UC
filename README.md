# Vedo2 Project

The repository for the code used in the Vedolizumab-Predict project and manuscript called "High-dimensional single-cell analysis identifies cellular signatures associated with response to vedolizumab therapy in ulcerative colitis".

The steps for replicating the results of this study are grouped into categories, and each in its own subdirectory including a README file.

# Overview 

For each different step of the scRNAseq analysis, separate codes were generated. Languages and packages are listed below: 

- R >= 3.6.1
- Seurat >= 3.1


External tools used were:

- Souporcell v.1: https://github.com/wheaton5/souporcell
- Cellranger v.3.0.2: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
- Azimuth v: https://azimuth.hubmapconsortium.org/
- fido v.1.0.2:
- ClusterProfiler v.4.4.4: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
- Cellchat v.1.4.0: http://www.cellchat.org/
- nichenetr v.1.1.0: https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md


the content for each step includes:

-   Alignment of the sequence data to the GRCh38 reference Human Genome
-   Demultiplexing, doublet_detection, and quality control (QC) 
-   General quality control and filtering
-   Normalization and data integration
-   Dimensional reduction and clustering
-   Cell type classification
-   Expression of integrins
-   Cell abundances analysis
-   Differential gene expression analysis and pathway analysis
-   Cell-cell interaction analysis

# Workflow

steps should be done in this order, to replicate the results

## Alignment of the sequence data to the GRCh38 reference Human Genome

The 1_Alignment contains the script 'Cellranger_script' which would create jobs to run each 10x lane through CellRanger.

## Demultiplexing, double detection, and quality control (QC)

