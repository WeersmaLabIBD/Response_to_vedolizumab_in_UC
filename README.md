# Vedolizumab-Predict

Repository of the code used in the Vedolizumab-Predict project.

The steps for replicating the results of this study are grouped into categories, and each in its own subdirectory including a README file.

# Contents

-   Alignment of the sequence data to the GRCh38 reference Human Genome
-   Demultiplexing_and_doublet_detection: demultiplexing and assignment of individuals to cells
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

## alignment

The 1_Alignment contains the script 'Cellranger script' which would create jobs to run each 10x lane through CellRanger.

## demultiplexing

6. '*demultiplexing_and_doublet_detection/lpmcv2_create_souping_samples_jobs.sh*' This will create Souporcell demultiplexing jobs. These jobs are then to be submitted. Be sure that there is an unzipped version of the barcodes in the counts directory

## merging

7. '*preprocessing/lpmcv2_lane_to_seurat.R*' This will take a 10x lane that is the output of cellranger, and turn it into a Seurat object

8. '*preprocessing/lpmcv2_create_lane_to_seurat_jobs.sh*' This script will create jobs that read a 10x lane and turn it into Seurat objects using the lane_to_seurat.R. These jobs then need to be submitted.

9. '*preprocessing/lpmcv2_filter_and_merge_lanes.R*' This will merge the Seurat objects for each lane, and combine them into one object. Then the Souporcell output is read, and doublets are removed. Finally the cells with a mitochondrial gene percentage of 70 percent or higher are removed

## normalization

10.  '*preprocessing/lpmcv2_sct_normalize_and_cluster.R*' perform SCT normalization, clustering and UMAP dimensional reduction

## cell type assignment

11. '*celltype_assigning/split_compartments.R*' This will do a compartment assignment of immune, stromal and epipithelial to the Seurat object, based on marker gene expression

12. '*celltype_assigning/lpmcv2_create_elmentaite2021_objects.R*' create a Seurat object from the Elmentaite 2021 study

13. '*celltype_assigning/lpmcv2_create_martin_object.R*' create a Seurat object from the Martin 2019 study

14. '*celltype_assigning/lpmcv2_add_azimuth_classifications.R' add the cell type assignments from Azimuth to the object

15. '*celltype_assigning/lpmcv2_add_lower_celltypes.R' add lower resolution cell types
