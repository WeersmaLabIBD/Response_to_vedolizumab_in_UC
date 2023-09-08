#############################################
# using azimuth to annotate vedo2 data      #
# change Smillie train data to fit our data #
#############################################

#############################################
# libraries                                 #
#############################################

library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(SeuratDisk)
library(readxl)
library(scPred)

#############################################
# functions                                 #
#############################################

# normalize, cluster and UMAP with the SCT data (for reference data)
normalize_and_cluster_sct_method_1 <- function(seurat_object, resolution=1.2, dims=1:30, k.param=20){
  # remove existing SCT normalization if present
  try({
    seurat_object@assays$SCT <- NULL
  })
  # start normalizing
  seurat_object <- SCTransform(seurat_object)
  # run PCA
  seurat_object <- RunPCA(seurat_object, verbose = TRUE)
  # calculate nearest neighbours
  seurat_object <- FindNeighbors(seurat_object, dims = dims, verbose = TRUE)
  # use same seed every time
  set.seed(0)
  # do clustering using knn
  seurat_object <- FindClusters(seurat_object, resolution = resolution, verbose = TRUE)
  # do 2d dimensional reduction
  seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = TRUE, return.model = TRUE)
  # return the result
  return(seurat_object)
}

# normalize, cluster and UMAP with the SCT data (for biopsy/query data)
normalize_and_cluster_sct_method_2 <- function(seurat_object, resolution=1.2, dims=1:30, k.param=20){
  # remove existing SCT normalization if present
  try({
    seurat_object@assays$SCT <- NULL
  })
  # start normalizing
  seurat_object <- SCTransform(seurat_object)
  # run PCA
  seurat_object <- RunPCA(seurat_object, verbose = TRUE)
  # calculate nearest neighbours
  seurat_object <- FindNeighbors(seurat_object, dims = dims, verbose = TRUE)
  # use same seed every time
  set.seed(0)
  # do clustering using knn
  seurat_object <- FindClusters(seurat_object, resolution = resolution, verbose = TRUE)
  # do 2d dimensional reduction
  seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = TRUE)
  # return the result
  return(seurat_object)
}

# find transfer anchors between the reference data and the query data, and map the query data to the reference data
azimuth_mapping <- function(seurat_object, reference){
  # find transfer anchors between the reference data and query data
  anchors <- FindTransferAnchors(
    reference = reference,
    query = seurat_object,
    normalization.method = "SCT",  
    reference.reduction = "pca",
    dims = 1:30,
    recompute.residuals = FALSE
  )
  # map the query data to the reference data
  seurat_object <- MapQuery(
    anchorset = anchors,
    query = seurat_object,
    reference = reference,
    refdata = list(
      cell_type.pred = "Cluster"
    ),
    reference.reduction = "pca",
    reduction.model = "umap"  
  )
  # return the result
  return(seurat_object)
}

#############################################
# main codes                                #
#############################################

# load the reference datasets
Smillie_epi_ref <- readRDS("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/cell-type-classification/machine-learning/reference/train.Epi.seur.rds")
Smillie_fib_ref <- readRDS("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/cell-type-classification/machine-learning/reference/train.Fib.seur.rds")
Smillie_imm_ref <- readRDS("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/cell-type-classification/machine-learning/reference/train.Imm.seur.rds")

# load our own dataset
vedo2 <- readRDS("/groups/umcg-weersma/tmp01/Emilia/batch_1/3_batch_1_integration/vedo2_batch1_biop_integrated_noribo_nonorm_withmito_ref.rds")

# read SingleR celltyping vedo metadata that Maria did
vedo2_biopsy_singleR_celltyping <- read_xlsx("/groups/umcg-weersma/tmp01/Frank/Azimuth_celltyping/singleR/Maria_Labels_AutomatedCellAnnotation.xlsx", sheet = 1, col_names = T)
vedo2_biopsy_singleR_celltyping <- as.data.frame(vedo2_biopsy_singleR_celltyping)
vedo2_biopsy_singleR_celltyping <- vedo2_biopsy_singleR_celltyping[, c("cell.names", "label_res_10")]

# change column name 
colnames(vedo2_biopsy_singleR_celltyping)[2] <- "singleR.annotation"

# add singleR cell annotation into vedo2 Seurat object
singleR.annotation <- vedo2_biopsy_singleR_celltyping$singleR.annotation
names(singleR.annotation) <- rownames(vedo2@meta.data)
vedo2 <- AddMetaData(
  object = vedo2,
  metadata = singleR.annotation,
  col.name = "singleR.annotation"
)

# add lineages (Epithelial, Immune, Stromal) as meta data  
# set defaultassay as "integrated"
DefaultAssay(vedo2) <- "integrated"

# perform linear dimensional reduction
vedo2 <- RunPCA(vedo2, verbose = FALSE)

# cluster the cells
set.seed(0)
vedo2 <- FindNeighbors(vedo2, dims = 1:30)
vedo2 <- FindClusters(vedo2, resolution = 0.7)

# run non-linear dimensional reduction (UMAP)  
vedo2 <- RunUMAP(vedo2, dims = 1:30)
DimPlot(vedo2, label = T) #32 clusters

# identify each comaprtment (epithelial cells, stromal cells, and immune cells) manually before
## epithelial cells (cluster 3, 5, 7, 12, 14, 16, 19, 27, 28)
## stromal cells (cluster 1, 9, 20, 25, 26, 29)
## immune cells (cluster 0, 2, 4, 6, 8, 10, 11, 13, 15, 17, 18, 21, 22, 23, 24, 30, 31)

# subset the dataset to the already classified coarse cell types
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==3),"lineage"]<-"Epithelial"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==5),"lineage"]<-"Epithelial"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==7),"lineage"]<-"Epithelial"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==12),"lineage"]<-"Epithelial"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==14),"lineage"]<-"Epithelial"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==16),"lineage"]<-"Epithelial"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==19),"lineage"]<-"Epithelial"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==27),"lineage"]<-"Epithelial"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==28),"lineage"]<-"Epithelial"

vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==1),"lineage"]<-"Stromal"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==9),"lineage"]<-"Stromal"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==20),"lineage"]<-"Stromal"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==25),"lineage"]<-"Stromal"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==26),"lineage"]<-"Stromal"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==29),"lineage"]<-"Stromal"

vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==0),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==2),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==4),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==6),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==8),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==10),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==11),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==13),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==15),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==17),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==18),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==21),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==22),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==23),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==24),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==30),"lineage"]<-"Immune"
vedo2@meta.data[which(vedo2@meta.data$seurat_clusters==31),"lineage"]<-"Immune"

vedo2_epi <- vedo2[, vedo2@meta.data[['lineage']] == 'Epithelial'] 
vedo2_stro <- vedo2[, vedo2@meta.data[['lineage']] == 'Stromal'] 
vedo2_imm <- vedo2[, vedo2@meta.data[['lineage']] == 'Immune'] 

# put the vedo objects through the pipeline
vedo2_epi <- normalize_and_cluster_sct_method_2(vedo2_epi)
vedo2_stro <- normalize_and_cluster_sct_method_2(vedo2_stro)
vedo2_imm <- normalize_and_cluster_sct_method_2(vedo2_imm)

# do the same normalization for reference data
Smillie_epi_ref <- normalize_and_cluster_sct_method_1(Smillie_epi_ref)
Smillie_stro_ref <- normalize_and_cluster_sct_method_1(Smillie_fib_ref)
Smillie_imm_ref <- normalize_and_cluster_sct_method_1(Smillie_imm_ref)

# make dimplots for Smillie epithelial refernece data
p1 <- DimPlot(smillie_epi_ref, group.by = "Cluster", label = TRUE, repel = TRUE) 
p2 <- DimPlot(smillie_epi_ref, group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
p1 + p2 
ggsave("Smillie_epi_train_cluster_celltypes.jpg", width= 20, height= 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

# make dimplots for Smillie stromal refernece data
p3 <- DimPlot(smillie_stro_ref, group.by = "Cluster", label = TRUE, repel = TRUE) 
p4 <- DimPlot(smillie_stro_ref, group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
p3 + p4 
ggsave("Smillie_stro_train_cluster_celltypes.jpg", width= 20, height= 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

# make dimplots for Smillie immune refernece data
p5 <- DimPlot(smillie_imm_ref, group.by = "Cluster", label = TRUE, repel = TRUE) 
p6 <- DimPlot(smillie_imm_ref, group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
p5 + p6 
ggsave("Smillie_imm_train_cluster_celltypes.jpg", width= 20, height= 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

# if the references have decent overlap of clusters versus annotation, we can do the azimuth stuff
# the answer is yes

################################################
# azimuth annotation on epithelial cells       #
################################################

vedo2_epi <- azimuth_mapping(seurat_object = vedo2_epi, reference = Smillie_epi_ref)
saveRDS(vedo2_epi, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/object/vedo2_batch1_epi_SCT_azimuth.rds")

# grab predicted metadata
vedo2_epi_metadata <- vedo2_epi@meta.data
vedo2_epi_metadata <- as.data.frame(vedo2_epi_metadata)
vedo2_epi_metadata <- vedo2_epi_metadata[, c("seurat_clusters", "singleR.annotation", "predicted.cell_type.pred.score", "predicted.cell_type.pred")]

# change column name 
colnames(vedo2_epi_metadata)[3] <- "azimuth prediction score"
colnames(vedo2_epi_metadata)[4] <- "azimuth annotation"
write.csv(vedo2_epi_metadata, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/metadata/vedo2_epi_metadata.csv")

# make dimplots for epithelial between azimuth and singleR
P1 <- DimPlot(vedo2_epi, group.by = "singleR.annotation", label = TRUE, repel = TRUE)
P2 <- DimPlot(vedo2_epi, group.by = "predicted.cell_type.pred", label = TRUE, repel = TRUE)
P1 + P2
ggsave("vedo2_epi_singleR_azimuth_celltypes.jpg", width= 20, height= 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

# compare cell number overlap of epithelial between azimuth and singleR
epi_overlap_cell_number_azimuth_singleR <- crossTab(vedo2_epi, "singleR.annotation", "predicted.cell_type.pred")
write.csv(epi_overlap_cell_number_azimuth_singleR, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/epi_overlap_cell_number_azimuth_singleR.csv")

# compare cell percentage overlap of epithelial between azimuth and singleR
epi_overlap_cell_percentage_azimuth_singleR <- crossTab(vedo2_epi, "singleR.annotation", "predicted.cell_type.pred", output = "prop")
write.csv(epi_overlap_cell_percentage_azimuth_singleR, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/epi_overlap_cell_percentage_azimuth_singleR.csv")

# make dimplots for epithelial between azimuth annotation and Seurat cluster
P3 <- DimPlot(vedo2_epi, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
P4 <- DimPlot(vedo2_epi, group.by = "predicted.cell_type.pred", label = TRUE, repel = TRUE)
P3 + P4
ggsave("vedo2_epi_clusters_azimuth.jpg", width= 20, height= 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

# compare overlap cell number between azimuth annotation and Seurat cluster
epi_overlap_cell_number_azimuth_cluster <- crossTab(vedo2_epi, "seurat_clusters", "predicted.cell_type.pred")
write.csv(epi_overlap_cell_number_azimuth_cluster, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/epi_overlap_cell_number_azimuth_cluster.csv")

# compare cell annotation between azimuth annotation and Seurat cluster
epi_overlap_cell_percentage_azimuth_cluster <- crossTab(vedo2_epi, "seurat_clusters", "predicted.cell_type.pred", output = "prop")
write.csv(epi_overlap_cell_percentage_azimuth_cluster, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/epi_overlap_cell_percentage_azimuth_cluster.csv")

# check the azimuth prediction score
# copy the prediction column in the metadata to a new column
vedo2_epi@meta.data$predicted.cell_type.pred.wreject <- as.character(vedo2_epi@meta.data$predicted.cell_type.pred)
# for all rows in metadata, where l1 score < 0.75, change the new column column to NA
vedo2_epi@meta.data[vedo2_epi@meta.data$predicted.cell_type.pred.score < 0.75, 'predicted.cell_type.pred.wreject'] <- NA
# plot the prediction cell type without changes
P5 <- DimPlot(vedo2_epi, reduction = 'umap', group.by = 'predicted.cell_type.pred', label = TRUE, label.box = T, repel = T)
# plot the prediction cell type, with the rejections (< 0.75 turned into NA)
P6 <- DimPlot(vedo2_epi, reduction = 'umap', group.by = 'predicted.cell_type.pred.wreject')
P5 + P6
# save the plot
ggsave('epi_predicted_clusters.jpg', width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

################################################
# azimuth annotation on stromal cells          #
################################################
vedo2_stro <- azimuth_mapping(seurat_object = vedo2_stro, reference = Smillie_stro_ref)
saveRDS(vedo2_stro, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/object/vedo2_batch1_stro_SCT_azimuth.rds")

# grab predicted metadata
vedo2_stro_metadata <- vedo2_stro@meta.data
vedo2_stro_metadata <- as.data.frame(vedo2_stro_metadata)
vedo2_stro_metadata <- vedo2_stro_metadata[, c("seurat_clusters", "singleR.annotation", "predicted.cell_type.pred.score", "predicted.cell_type.pred")]

# change column name 
colnames(vedo2_stro_metadata)[3] <- "azimuth prediction score"
colnames(vedo2_stro_metadata)[4] <- "azimuth annotation"
write.csv(vedo2_stro_metadata, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/metadata/vedo2_stro_metadata.csv")

# make dimplots for epithelial between azimuth and singleR
P7 <- DimPlot(vedo2_stro, group.by = "singleR.annotation", label = TRUE, repel = TRUE)
P8 <- DimPlot(vedo2_stro, group.by = "predicted.cell_type.pred", label = TRUE, repel = TRUE)
P7 + P8
ggsave("vedo2_stro_singleR_azimuth_celltypes.jpg", width= 20, height= 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

#compare cell number overlap of stromal between azimuth and singleR
stro_overlap_cell_number_azimuth_singleR <- crossTab(vedo2_stro, "singleR.annotation", "predicted.cell_type.pred")
write.csv(stro_overlap_cell_number_azimuth_singleR, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/stro_overlap_cell_number_azimuth_singleR.csv")

#compare cell percentage overlap of stromal between azimuth and singleR
stro_overlap_cell_percentage_azimuth_singleR <- crossTab(vedo2_stro, "singleR.annotation", "predicted.cell_type.pred", output = "prop")
write.csv(stro_overlap_cell_percentage_azimuth_singleR, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/stro_overlap_cell_percentage_azimuth_singleR.csv")

# make dimplots for stromal between azimuth annotation and Seurat cluster
P9 <- DimPlot(vedo2_stro, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
P10 <- DimPlot(vedo2_stro, group.by = "predicted.cell_type.pred", label = TRUE, repel = TRUE)
P9 + P10
ggsave("vedo2_stro_clusters_azimuth.jpg", width= 20, height= 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

#compare overlap cell number between azimuth and Seurat cluster
stro_overlap_cell_number_azimuth_cluster <- crossTab(vedo2_stro, "seurat_clusters", "predicted.cell_type.pred")
write.csv(stro_overlap_cell_number_azimuth_cluster, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/stro_overlap_cell_number_azimuth_cluster.csv")

#compare cell annotation between azimuth and Seurat cluster
stro_overlap_cell_percentage_azimuth_cluster <- crossTab(vedo2_stro, "seurat_clusters", "predicted.cell_type.pred", output = "prop")
write.csv(stro_overlap_cell_percentage_azimuth_cluster, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/stro_overlap_cell_percentage_azimuth_cluster.csv")

# check the azimuth prediction score
# copy the prediction column in the metadata to a new column
vedo2_stro@meta.data$predicted.cell_type.pred.wreject <- as.character(vedo2_stro@meta.data$predicted.cell_type.pred)
# for all rows in metadata, where l1 score < 0.75, change the new column column to NA
vedo2_stro@meta.data[vedo2_stro@meta.data$predicted.cell_type.pred.score < 0.75, 'predicted.cell_type.pred.wreject'] <- NA
# plot the prediction cell type without changes
P11 <- DimPlot(vedo2_stro, reduction = 'umap', group.by = 'predicted.cell_type.pred', label = TRUE, label.box = T, repel = T)
# plot the prediction cell type, with the rejections (< 0.75 turned into NA)
P12 <- DimPlot(vedo2_stro, reduction = 'umap', group.by = 'predicted.cell_type.pred.wreject')
P11 + P12
# save the plot
ggsave('stro_predicted_clusters.jpg', width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

################################################
# azimuth annotation on immune cells           #
################################################

vedo2_imm <- azimuth_mapping(seurat_object = vedo2_imm, reference = Smillie_imm_ref)
saveRDS(vedo2_imm, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/object/vedo2_batch1_imm_SCT_azimuth.rds")

# grab predicted metadata
vedo2_imm_metadata <- vedo2_imm@meta.data
vedo2_imm_metadata <- as.data.frame(vedo2_imm_metadata)
vedo2_imm_metadata <- vedo2_imm_metadata[, c("seurat_clusters", "singleR.annotation", "predicted.cell_type.pred.score", "predicted.cell_type.pred")]

# change column name 
colnames(vedo2_imm_metadata)[3] <- "azimuth prediction score"
colnames(vedo2_imm_metadata)[4] <- "azimuth annotation"
write.csv(vedo2_imm_metadata, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/metadata/vedo2_imm_metadata.csv")

# make dimplots for epithelial between azimuth and singleR
P13 <- DimPlot(vedo2_imm, group.by = "singleR.annotation", label = TRUE, repel = TRUE)
P14 <- DimPlot(vedo2_imm, group.by = "predicted.cell_type.pred", label = TRUE, repel = TRUE)
P13 + P14
ggsave("vedo2_imm_singleR_azimuth_celltypes.jpg", width= 20, height= 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

#compare cell number overlap of immune between azimuth and singleR
imm_overlap_cell_number_azimuth_singleR <- crossTab(vedo2_imm, "singleR.annotation", "predicted.cell_type.pred")
write.csv(imm_overlap_cell_number_azimuth_singleR, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/imm_overlap_cell_number_azimuth_singleR.csv")

#compare overlap cell percentage overlap of immune between azimuth and singleR
imm_overlap_cell_percentage_azimuth_singleR <- crossTab(vedo2_imm, "singleR.annotation", "predicted.cell_type.pred", output = "prop")
write.csv(imm_overlap_cell_percentage_azimuth_singleR, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/imm_overlap_cell_percentage_azimuth_singleR.csv")

# make dimplots for stromal between azimuth annotation and Seurat cluster
P15 <- DimPlot(vedo2_imm, group.by = "seurat_clusters", label = TRUE, repel = TRUE)
P16 <- DimPlot(vedo2_imm, group.by = "predicted.cell_type.pred", label = TRUE, repel = TRUE)
P15 + P16
ggsave("vedo2_imm_clusters_azimuth.jpg", width= 20, height= 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")

#compare overlap cell number between azimuth and Seurat cluster
imm_overlap_cell_number_azimuth_cluster <- crossTab(vedo2_imm, "seurat_clusters", "predicted.cell_type.pred")
write.csv(imm_overlap_cell_number_azimuth_cluster, "//groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/imm_overlap_cell_number_azimuth_cluster.csv")

#compare cell annotation between azimuth and Seurat cluster
imm_overlap_cell_percentage_azimuth_cluster <- crossTab(vedo2_imm, "seurat_clusters", "predicted.cell_type.pred", output = "prop")
write.csv(imm_overlap_cell_percentage_azimuth_cluster, "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/overlap_singleR_azimuth/stro_overlap_cell_percentage_azimuth_cluster.csv")

# check the azimuth prediction score
# copy the prediction column in the metadata to a new column
vedo2_imm@meta.data$predicted.cell_type.pred.wreject <- as.character(vedo2_imm@meta.data$predicted.cell_type.pred)
# for all rows in metadata, where l1 score < 0.75, change the new column column to NA
vedo2_imm@meta.data[vedo2_imm@meta.data$predicted.cell_type.pred.score < 0.75, 'predicted.cell_type.pred.wreject'] <- NA
# plot the prediction cell type without changes
P17 <- DimPlot(vedo2_imm, reduction = 'umap', group.by = 'predicted.cell_type.pred', label = TRUE, label.box = T, repel = T)
# plot the prediction cell type, with the rejections (< 0.75 turned into NA)
P18 <- DimPlot(vedo2_imm, reduction = 'umap', group.by = 'predicted.cell_type.pred.wreject')
P17 + P18
# save the plot
ggsave('imm_predicted_clusters.jpg', width = 20, height = 10, path = "/groups/umcg-weersma/tmp01/Frank/batch1_vedo2_biopsy/Azimuth_celltyping/SCT_normalization/plots/")