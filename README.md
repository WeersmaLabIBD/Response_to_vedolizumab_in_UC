# Vedo2 Project

The repository for the code used in the Vedolizumab-Predict project and manuscript called "High-dimensional single-cell analysis identifies cellular signatures associated with response to vedolizumab therapy in ulcerative colitis".

The steps for replicating the results of this study are grouped into categories, and each in its own subdirectory including a README file.

# Overview 

For each different step of the scRNAseq analysis, separate codes were generated. Languages and packages are listed below: 

- R >= 3.6.1
- Seurat >= 3.1
- Python >= 3.7.4
- numpy 1.19.5
- pandas 1.2.1


External tools used were:

- Souporcell v.1: https://github.com/wheaton5/souporcell
- Cellranger v.3.0.2: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
- fido v.1.0.2: https://jsilve24.github.io/fido/articles/introduction-to-fido.html
- scCODA v.0.1.9 https://github.com/theislab/scCODA/blob/master/docs/source/getting_started.ipynb
- ClusterProfiler v.4.4.4: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
- Cellchat v.1.4.0: http://www.cellchat.org/
- nichenetr v.1.1.0: https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md

# Workflow
the content for each step includes these steps in this order:

-   Alignment of the sequence data to the GRCh38 reference Human Genome: https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/1_Alignment
-   Lane integration: https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/2_Lane_integration
-   Demultiplexing, doublet_detection, and quality control (QC): https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/3_Demultiplexing_and_doublet_detection
-   General quality control and filtering: https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/4_Preprocessing_and_celltype_annotation
-   Normalization and data integration: https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/4_Preprocessing_and_celltype_annotation
-   Dimensional reduction and clustering: https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/4_Preprocessing_and_celltype_annotation
-   Cell type classification: https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/4_Preprocessing_and_celltype_annotation
-   Expression of integrins: https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/5_integrin_expression
-   Cell abundances analysis: https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/6_Cell_abundances_analysis
-   Differential gene expression analysis and pathway analysis: https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/7_Differential_expression_analysis_and_Pathway
-   Cell-cell interaction analysis: https://github.com/WeersmaLabIBD/Response_to_vedolizumab_in_UC/tree/main/8_Cell2Cell_communication_analysis


# Hardware

Analyses were performed on either a 2021 MacBook Pro (32GB), the Gearshift cluster http://docs.gcc.rug.nl/gearshift/cluster/ 

