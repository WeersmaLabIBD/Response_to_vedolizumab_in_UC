# Vedolizumab-Predict

Repository of the code used in the Vedolizumab-Predict project.

The steps for replicating the results of this study are grouped into categories, and each in its own subdirectory including a README file.

# Contents

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

