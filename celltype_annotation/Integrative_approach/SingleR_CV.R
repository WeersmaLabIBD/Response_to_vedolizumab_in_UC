library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(dplyr)
library(data.table)
library(stats)
library(scran)
library(SingleR)

source("functions_util.R")

test_objects = readRDS("files/Smillie_Stromal_folds/test_objects.RDS")
train_objects = readRDS("files/Smillie_Stromal_folds/train_objects.RDS")



## Running SingleR algorithm. It has three versions:
# 1 - without markers (default)
# 2 - with list of markers corresponding to each cell type (add markers)
# 3 - with DE markers corresponding to each pairwise comparisons (add DEGs)


# create a list of cell-specific markers imported directly from Smillie
markers_smillie_all = read.delim("files/Fibro_markers.csv",header = T, sep = ",")

smillie_all=colnames(markers_smillie_all)
smillie_markers = as.list(smillie_all)
names(smillie_markers) = smillie_all

for (marker in names(smillie_markers)){
  mgenes= markers_smillie_all[[marker]][markers_smillie_all[[marker]]!=""]
  # check if these markers are exist in query dataset
  # mgenes = mgenes[mgenes %in% rownames(query_sce)]
  smillie_markers[[marker]] = mgenes
}

names(smillie_markers) =  c("Myofibroblasts", "WNT2B+ Fos-hi",  "WNT2B+ Fos-lo 1",  "WNT2B+ Fos-lo 2","RSPO3+","WNT5B+ 1","WNT5B+ 2", "Inflammatory Fibroblasts","Endothelial","Microvascular","Post-capillary Venules", "Pericytes")



# create a list of pruned markers for including DEGs into SingleR annotation
 geneSets_flt_markers <- fread("files/geneSets_markers_Smillie_Stroma.txt")


# Running each of these SingleR methods or all sequentially //Better to be done 3 times in-parallel or stick to the DEG version 
default=T
add_markers = F
add_DEGs = F
ref_assay ="RNA"
query_assay = "RNA"



#### assigning input parameters properly based on the version of SingleR which needs to be run (you can run one of the SingleR versions by setting others to False):
# 1 - Default version of SingleR : Set to False (default=False) if you want to skip this run
# 2 - SingleR with input list of celltype-specific markers per celltype: Set to False (add_markers=False) if you want to skip this version, otherwise make sure to assign "marker_list" based on the example file 
# 3 - SingleDE with input list of DE markers: Setto False (add_DEGs= False) if you want to skip this version, otherwise make sure to assign "DEGs_list" based on the example file
# please make sure at least one of the SingleR versions is properly set to be run...

# Cluster: is a name of the column in reference Seurat object which includes the reference cell labels


for (i in c(1:length(test_objects))){
  print(paste0("Running SingleR on object ",i))
  test_objects[[i]] = Run_SingleR(train_objects[[i]], test_objects[[i]],ref_assay,
                          query_assay,   default ,add_markers,
                          smillie_markers, add_DEGs, 
                          geneSets_flt_markers, "Cluster")
  
}


saveRDS(test_objects,"files/Smillie_Stromal_folds/test_objects.RDS")
