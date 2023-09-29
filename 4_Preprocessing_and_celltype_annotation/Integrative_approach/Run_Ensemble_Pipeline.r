# ## R Script
# ---
# Title: "Ensemble framework for scRNA-seq annotation"
# Author: Author: [Saman Farahmand] (mailto:saman.farahmand@takeda.com) - GI-DDU
# Date: "02/24/2022"
# Output: Seurat object in RDS format


# ---

# Description

#

library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(scmap)
library(SingleR)
library(dplyr)
library(scclusteval)
library(SCINA)
library(data.table)
library(scran)
library(caret)
library(singleCellNet)
source("/mnt/efs/prj/Saman-Maria/Groningen/Vedopredict/mfunctions.R")

flag_preproc = T
flag_SCINA = T
# Celltype-specific markers must be loaded if [flag_SCINA = TRUE] - Use the sample file to create the table
# Make sure the labels of celltypes are compatible with the cell-type labels in the reference data (Cluster column in the metadata)
SCINA_markers = read.delim("~/Projects/Groningen/Vedopredict/files/immune_markers.txt",header = T, sep = "\t")
flag_scmap_cell = F
flag_scmap_cluster = F

#### assigning input parameters properly based on the version of SingleR which needs to be run (you can run one of the SingleR versions by setting others to False):
# 1 - Default version of SingleR : Set to False (default=False) if you want to skip this run
# 2 - SingleR with input list of celltype-specific markers per celltype: Set to False (add_markers=False) if you want to skip this version, otherwise make sure to assign "marker_list" based on the example file 
# 3 - SingleDE with input list of DE markers: Set to False (add_DEGs= False) if you want to skip this version, otherwise make sure to assign "DEGs_list" based on the example file
# please make sure at least one of the SingleR versions is properly set to be run...


flag_SR_default=F
flag_SR_add_markers = F
# create a list of cell-specific markers imported directly from Smillie
# Celltype-specific markers must be loaded if [flag_SCINA = TRUE] - Use the sample file to create the table
# Make sure the labels of celltypes are compatible with the cell-type labels in the reference data (Cluster column in the metadata)
SR_celltype_markers = read.delim("~/Projects/Groningen/Vedopredict/files/immune_markers.txt",header = T, sep = "\t")


flag_SR_add_DEGs = T
# create a list of pruned markers for including DEGs into SingleR annotation
# Celltype-specific markers must be loaded if [flag_SCINA = TRUE] - Use the sample file to create the table
# Make sure the labels of celltypes are compatible with the cell-type labels in the reference data (Cluster column in the metadata)
#cellTypes_2_exclude = c("Glia")
 cellTypes_2_exclude = ""
geneSets_flt_markers <- fread("~/Projects/Groningen/Vedopredict/files/geneSets_markers_Smillie_Immune.txt")
#geneSets_flt_markers <- fread("~/Projects/Groningen/Vedopredict/files/geneSets_markers_Smillie_Stroma_old.txt")


# flag to set Azimuth method for running
flag_azimuth = TRUE

# flag to set SCT method for Running
flag_SCN = TRUE

# Selected methods for consensus annotation
# - SCINA
# - scmap_cluster
# - scmap_cell
# - singleR [for the default version]
# - singleRm [for the version with input marker list]
# - singleRDE [for the version with DE markers]
# - azimuth
# - SCN
methods_cons = c("singleRDE","azimuth","SCN")


# Defining output directory and output RDS file names
output_dir = "~/Projects/Groningen/Vedopredict/Vedolizumab/"
filename = "Vedo_Immune_annot.RDS"


### Loading the reference and query Seurat objects - This script accepts inputs as RDS formatted files including Seurat objects

# It's assumed that the objects are not normalized and scaled. It's recommended to apply same normalization, scaling and clustering approaches on both objects for consistency
# In case you have already normalized your data, you can skip this step by setting [flag_preproc = FALSE]. We use SCTransform for normalization.
# Make sure the label column for cell types in reference Seurat metadata is named [Cluster]

ref_seur = readRDS("~/Projects/Groningen/Vedopredict/Smilie/github/train.Imm.seur.rds")
quer_seur = readRDS("~/Projects/Groningen/Vedopredict/Vedolizumab/Vedo_immune.RDS")

# normalizing and clustering the reference and query objects
if(flag_preproc) {
  
  print("Normalizing and clustering the reference dataset...")
  ref_seur = normalize_and_cluster_sct_method(ref_seur, resolution = 0.8, dims = 1:30) 
  print("Normalizing and clustering the query dataset...")
  quer_seur = normalize_and_cluster_sct_method(quer_seur, resolution = 0.8, dims = 1:30)
}

### Running SCINA method

if(flag_SCINA) {
  print("Running SCINA method...")
  ## Format the list of markers
  label_all=SCINA_markers$cell_type
  markers = as.list(label_all)
  names(markers) = label_all
  for (marker in label_all){
    markers[[marker]] = unlist(strsplit(SCINA_markers[SCINA_markers$cell_type==marker,]$markers, split = ","))
  }
  # Running SCINA prediction
  quer_seur = Run_SCINA(quer_seur, markers = markers)

  }


### Running SCmap method

if(flag_scmap_cell | flag_scmap_cluster) {
  
  quer_seur = Run_Scmap(ref_seur,quer_seur,ref_assay ="RNA", query_assay = "RNA",
                                scmap_cell = flag_scmap_cell,
                                scmap_cluster = flag_scmap_cluster,
                                n_features=500, Cluster="Cluster")
  
}


### Running SingleR method

if(flag_SR_default | flag_SR_add_markers | flag_SR_add_DEGs) {
  
  print("Running SingleR method...")
  quer_seur = Run_SingleR(ref_seur, quer_seur,ref_assay ="RNA",
                                    query_assay = "RNA", flag_SR_default ,flag_SR_add_markers,
                                    SR_celltype_markers, flag_SR_add_DEGs, cellTypes_2_exclude,
                                    geneSets_flt_markers, Cluster = "Cluster")
  
}

### Running Azimuth method

if(flag_azimuth) {
  
  
  print("Running Azimuth method...")
  ref_seur = SCT_umap (ref_seur, resolution=1.2, dims=1:30)
  quer_seur = SCT_umap (quer_seur , resolution=1.2, dims=1:30)
  #The reference was normalized using SCT, so we use the same approach to normalize the query here.
  re.obj = azimuth_mapping(quer_seur, ref_seur,"Cluster","SCT")
  
  quer_seur$azimuth =  re.obj$predicted.cell_type.pred
  quer_seur$Azimuth_scores = re.obj$predicted.cell_type.pred.score
  
  
  
  
}


### Running SCN method
if (flag_SCN) {
  
  print("Running SCN method...")
  # convert query object to SCN object
  # exp_type options can be: counts, normcounts, and logcounts, if they are available in your object. Default is counts
  quer_seur$cellid = rownames(quer_seur@meta.data)
  seur_query = extractSeurat(quer_seur, exp_slot_name = "counts")
  quer.sampTab = seur_query$sampTab
  quer.expDat = seur_query$expDat
  
  # convert reference object to SCN object
  
  seur_ref = extractSeurat(ref_seur, exp_slot_name = "counts")
  ref.sampTab = seur_ref$sampTab
  ref.expDat = seur_ref$expDat
  ref.sampTab$cell = rownames(ref.sampTab)
  
  commonGenes = intersect(rownames(ref.expDat), rownames(quer.expDat))

  ref.expDat = ref.expDat[commonGenes,]
  quer.expDat = quer.expDat[commonGenes,]
  
  
  print("Training SCN model using the reference data...")
  class_info = scn_train(stTrain = ref.sampTab, expTrain = ref.expDat, nTopGenes = 30, nRand = 70,
                         nTrees = 1000, nTopGenePairs = 60, dLevel = "Cluster", colName_samp = "cell")
  
  nqRand = 50
  quer_predict = scn_predict(class_info[['cnProc']], quer.expDat, nrand=nqRand)
  
  
  # This classifies a cell with  the catgory with the highest classification score or higher than a classification score threshold of your choosing.
  # The annotation result can be found in a column named category in the query sample table.
  
  stPark <- get_cate(classRes = quer_predict, sampTab = quer.sampTab, dLevel = "azimuth", sid = "cellid", nrand = nqRand)
  
  quer_seur$SCN = stPark$category
  quer_seur$SCN_scores = stPark$scn_score
  quer_seur$cellid = NULL
}

print("Performing the consensus annotation for ")
quer_seur$consensus = apply(quer_seur@meta.data[,methods_cons], 1, get_consensus_label)

saveRDS(quer_seur,paste0(output_dir,filename))


### Running SCN method