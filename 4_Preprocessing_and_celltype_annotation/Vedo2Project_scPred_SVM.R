########################################
# using scPred to annotate  vedo data #
########################################

########################################
# libraries                            #
########################################

library("scPred")
library("Seurat")
library("magrittr")
library("ggplot2")
library("lattice")
library("readxl")

########################################
# functions                            #
########################################

### Normalize, cluster and UMAP with the SCT data
normalize_and_cluster_sct_method <- function(seurat_object, resolution=1.2, dims=1:30, k.param=20){
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
  set.seed(1234)
  # do clustering using knn
  seurat_object <- FindClusters(seurat_object, resolution = resolution, verbose = TRUE)
  # do 2d dimensional reduction
  seurat_object <- RunUMAP(seurat_object, dims = 1:30, verbose = TRUE)
  # return the result
  return(seurat_object)
}

### Change, probabilities

setMethod("get_probabilities", 
          signature("scPred"), 
          function(object){
            
            models <- get_classifiers(object)
            
            probs <- lapply(models, function(x){
              i <- x$levels[x$levels != "other"]
              x$pred[c(i, "rowIndex")]
              
            })
            
            barcodes <- rownames(get_metadata(object))
            #print(length(barcodes))
            
            
            probs <- mapply(function(x, x_name){
              res <- data.frame(x, barcodes = barcodes[x$rowIndex])
              res$rowIndex <- NULL
              names(res)[1] <- x_name
              res
            }, probs, names(probs), SIMPLIFY = FALSE)
            
            probs <- Reduce(function(x, y) merge(x, y, by = "barcodes"), probs)
            
            #print(probs)
            #print(nrow(probs))
            
            # strange error that we are hacking around
            if(nrow(probs) != length(unique(probs$barcodes))){
              # though we should at least warn them
              warning(paste('prob row numbers is', as.character(nrow(probs)),
                            ', while unique barcode number is', as.character(length(unique(probs$barcodes))),
                            ', with original metadata barcode being', as.character(length(barcodes)),
                            ', and barcode uniqueness being', as.character(length(unique(barcodes))),
                            ', causing duplicates to be removed'
              )
              )
              # get the cells that are duplicated
              duplicates <- probs[duplicated(probs$barcodes), 'barcodes']
              # print the duplicates
              print(duplicates)
              probs <- probs[match(as.character(barcodes), probs$barcodes), ]
            }
            
            rownames(probs) <- probs$barcodes
            probs$barcodes <- NULL
            
            probs <- probs[match(barcodes, rownames(probs)), ]
            probs
            
          })

########################################
# Main code Normalization              #
########################################

#Location of references
reference_Smillie_Epi_loc <- ("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/cell-type-classification/machine-learning/reference/train.Epi.seur.rds")
reference_Smillie_Stro_loc <- ("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/cell-type-classification/machine-learning/reference/train.Fib.seur.rds")
reference_Smillie_Imm_loc <- ("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/cell-type-classification/machine-learning/reference/train.Imm.seur.rds")
#Location of the query dataset
object_loc_all <- ("/groups/umcg-weersma/tmp01/Emilia/batch_1/3_batch_1_integration/vedo2_batch1_biop_integrated_noribo_nonorm_withmito_ref.rds")

#Load the references
smillie_epi_ref <- readRDS(reference_Smillie_Epi_loc)
smillie_stro_ref <- readRDS(reference_Smillie_Stro_loc)
smillie_imm_ref <- readRDS(reference_Smillie_Imm_loc)

#Load the query dataset
vedo2 <- readRDS(object_loc_all)

#Read Single R celltyping vedo metadata 
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

#Add lineages (Epithelial, Immune, Stromal) as meta data based on resolution 0.7 to get 32 clusters
vedo2 = UpdateSeuratObject(object = vedo2)
DefaultAssay(vedo2) <- "integrated"
set.seed(1234)
vedo2 <- FindClusters(vedo2, resolution = 0.7)
DimPlot(vedo2, label = T) 

#32 clusters
# Epi	3	5	7	12	14	16	19	27	28								
# Stromal	1	9		20	25	26	29										
# Immune	0	2	4	6	8	10	11	13 15	17	18	21	22	23	24	30	31

#Add metadata of the three compartments (manual cluster based annotation)           
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

# subset the dataset to the already classified coarse cell types
vedo2_epi <- vedo2[, vedo2@meta.data[['lineage']] == 'Epithelial'] 
vedo2_stro <- vedo2[, vedo2@meta.data[['lineage']] == 'Stromal'] 
vedo2_imm <- vedo2[, vedo2@meta.data[['lineage']] == 'Immune'] 

#Normalization
vedo2_epi <- normalize_and_cluster_sct_method(vedo2_epi)
vedo2_stro <- normalize_and_cluster_sct_method(vedo2_stro)
vedo2_imm <- normalize_and_cluster_sct_method(vedo2_imm)
# Save

#Normalization
smillie_epi_ref <- normalize_and_cluster_sct_method(smillie_epi_ref)
smillie_stro_ref <- normalize_and_cluster_sct_method(smillie_stro_ref)
smillie_imm_ref <- normalize_and_cluster_sct_method(smillie_imm_ref)
# Save

#Check epithelial cells
p1 <- DimPlot(smillie_epi_ref, group.by = "Cluster", label = TRUE, repel = TRUE) 
p2 <- DimPlot(smillie_epi_ref, group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
p1 + p2 

#Check stromal cells 
p1 <- DimPlot(smillie_stro_ref, group.by = "Cluster", label = TRUE, repel = TRUE) 
p2 <- DimPlot(smillie_stro_ref, group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
p1 + p2 

#Check immune cells
p1 <- DimPlot(smillie_imm_ref, group.by = "Cluster", label = TRUE, repel = TRUE) 
p2 <- DimPlot(smillie_imm_ref, group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
p1 + p2 

# if the references have decent overlap of clusters versus annotation, you can continue the scpred script



########################################
# scPred  on Epithelial cells          #
########################################

# Upload files 
smillie_epi_ref <- ("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/seurat_preprocess_samples/objects/NKscPred/smillie_epi_ref_scpred.rds")
vedo2_epi <- ("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/seurat_preprocess_samples/objects/NKscPred/vedo2_epi_scpred.rds")
smillie_epi_ref <- readRDS(smillie_epi_ref)
vedo2_epi <- readRDS(vedo2_epi)

#Training classifiers with scPred
smillie_epi_ref <- getFeatureSpace(smillie_epi_ref, "Cluster") 
smillie_epi_ref <- trainModel(smillie_epi_ref)
# Save
# Run probabilities function script
get_probabilities(smillie_epi_ref) %>% head()
get_scpred(smillie_epi_ref)
plot_probabilities(smillie_epi_ref)
# Looks well trained

#Cell classification
vedo2_epi <- scPredict(vedo2_epi, smillie_epi_ref_trained
# Save
# Plot the classifications over the aligned data.
DimPlot(vedo2_epi, group.by = "scpred_prediction", reduction = "scpred")
# Cun UMAP using the aligned data as an input
vedo2_epi <- RunUMAP(vedo2_epi, reduction = "scpred", dims = 1:30) 
DimPlot(vedo2_epi, group.by = "scpred_prediction", reduction ="umap", label = TRUE, repel = TRUE) 
# Compare the results with the original labels
DimPlot(vedo2_epi, group.by = "singleR.annotation", reduction ="umap", label = TRUE, repel = TRUE) 

# Verify the performance of the models in the query dataset by using crossTab to create a contingency table of two colums from the metadata.
table_ann_epi <- crossTab(vedo2_epi, "singleR.annotation", "scpred_prediction")
table_ann_epi <- data.frame(table_ann_epi)

# The proportion of cells can be obtained using output = "prop
table_ann_epi_2 <- crossTab(vedo2_epi, "singleR.annotation", "scpred_prediction", output = "prop")
table_ann_epi_2 <- data.frame(table_ann_epi_2)



########################################
# scPred  on Stromal cells             #
########################################

# Upload files 
smillie_stro_ref <- ("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/seurat_preprocess_samples/objects/NKscPred/smillie_stro_ref_scpred.rds")
vedo2_stro <- ("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/seurat_preprocess_samples/objects/NKscPred/vedo2_stro_scpred.rds")
smillie_stro_ref <- readRDS(smillie_stro_ref)
vedo2_stro <- readRDS(vedo2_stro)

#Training classifiers with scPred
smillie_stro_ref <- getFeatureSpace(smillie_stro_ref, "Cluster")
smillie_stro_ref <- trainModel(smillie_stro_ref)
# Save
# Run probabilities function script
get_probabilities(smillie_stro_ref) %>% head()
get_scpred(smillie_stro_ref)
plot_probabilities(smillie_stro_ref)
# Looks well trained

#Cell classification
vedo2_stro <- scPredict(vedo2_stro, smillie_stro_ref_trained)
# Save
# Plot the classifications over the aligned data.
DimPlot(vedo2_stro, group.by = "scpred_prediction", reduction = "scpred")
# Run UMAP using the aligned data as an input
vedo2_stro <- RunUMAP(vedo2_stro, reduction = "scpred", dims = 1:30) #original labels
DimPlot(vedo2_stro, group.by = "scpred_prediction", label = TRUE, repel = TRUE) #predicted labels
# Compare the results with the original labels
DimPlot(vedo2_stro, group.by = "singleR.annotation", label = TRUE, repel = TRUE) #compare results

# Verify the performance of the models in the query dataset by using crossTab to create a contingency table of two colums from the metadata.
table_ann_stro<- crossTab(vedo2_stro, "singleR.annotation", "scpred_prediction")
table_ann_stro <- data.frame(table_ann_stro)

# The proportion of cells can be obtained using output = "prop
table_ann_stro_2 <- crossTab(vedo2_stro, "singleR.annotation", "scpred_prediction", output = "prop")
table_ann_stro_2 <- data.frame(table_ann_stro_2)


########################################
# scPred  on Immune cells              #
########################################
# Upload files again if necessary 
smillie_imm_ref <- ("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/seurat_preprocess_samples/objects/NKscPred/smillie_imm_ref_scpred.rds") 
vedo2_imm <- ("/groups/umcg-weersma/tmp01/projects/vedo2_predict/ongoing/seurat_preprocess_samples/objects/NKscPred/vedo2_imm_scpred.rds")
smillie_imm_ref <- readRDS(smillie_imm_ref)
vedo2_imm <- readRDS(vedo2_imm)

#Training classifiers with scPred
smillie_imm_ref <- getFeatureSpace(smillie_imm_ref, "Cluster")
smillie_imm_ref <- trainModel(smillie_imm_ref)
# Save
# Run probabilities function script
get_probabilities(smillie_imm_ref) %>% head()
get_scpred(smillie_imm_ref)
plot_probabilities(smillie_imm_ref)
# Looks well trained

#Cell classification
vedo2_imm <- scPredict(vedo2_imm, smillie_imm_ref)
# Save
# Plot the classifications over the aligned data.
DimPlot(vedo2_imm, group.by = "scpred_prediction", reduction = "scpred")
# Error scpred not found
# We cannot continue this part of the script.

# run UMAP using the aligned data as an input
vedo2_imm <- RunUMAP(vedo2_imm, reduction = "scpred", dims = 1:30) #original labels
DimPlot(vedo2_imm, group.by = "scpred_prediction", label = TRUE, repel = TRUE) #predicted labels
# compare the results with the original labels
DimPlot(vedo2_imm, group.by = "singleR.annotation", label = TRUE, repel = TRUE) #compare results

#Verify the performance of the models in the query dataset by using crossTab to create a contingency table of two colums from the metadata.
crossTab(vedo2_imm, "cell_type", "scpred_prediction")
table_ann_imm <- data.frame(table_ann_imm)

#The proportion of cells can be obtained using output = "prop
table_ann_imm <- crossTab(vedo2_imm, "cell_type", "scpred_prediction", output = "prop")
table_ann_imm <- data.frame(table_ann_imm)




