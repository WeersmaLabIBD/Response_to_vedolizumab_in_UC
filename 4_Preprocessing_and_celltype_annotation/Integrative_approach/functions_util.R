library(dplyr)
library(ComplexHeatmap)
library(circlize)

# clusters_col: name of a column including the cluster numbers 
# label_col column: name of a column including the cell annotations
# cellid column: name of a column of cell identifiers

# Return an annotated seurat object
Run_SCINA = function(seur_obj,markers){
  
  # Convert the expression data from Seurat object into a matrix data structure
  exprMatrix <- as.matrix(Seurat::GetAssayData(seur_obj))
  # Run SCINA on the query data using the marker genes to identify cell types
  # Specifying rm_overlap = FALSE allows the same marker gene to specify multiple cell types which may be useful if identifying cell subtypes or other similar types of cells
  # Specifying allow_unknown = TRUE allows cells to be labeled as "unknown" instead of being assigned a low-confident label
  predictions.scina = SCINA::SCINA(exp = exprMatrix, signatures = markers,
                                   rm_overlap = FALSE, allow_unknown = TRUE)
  # Add SCINA annotation information to each cell in Seurat object
  seur_obj$SCINA <- predictions.scina$cell_labels
  return(seur_obj)
}


## Run scmap approach to annotate Vedo object using Smillie Reference
# Scmap uses data formatted as a 'SingleCellExperiment object', and assumes by default that 
# gene names are found in a column named 'feature_symbol' while the cell-type labels are in a 
# column named 'cell_type1'. In addition, scmap requires that you normalize and log-transform the reference data
Run_Scmap = function(ref_seur_obj, query_seur_obj, ref_assay ="RNA", 
                     query_assay = "RNA", scmap_cell = T,
                     scmap_cluster = T, 
                     n_features, Cluster) {
  
  
  # ref_seur_obj = test_objects[[1]]
  # query_seur_obj = train_objects[[1]]
  # ref_assay ="RNA"
  # query_assay = "RNA"
  # n_features = 1000
  # Cluster = "Cluster"
  
  
  ref_sce = as.SingleCellExperiment(ref_seur_obj,ref_assay)
  colnames(colData(ref_sce))
  
  # Assign cell-type labels in a column named "cell_type1"
  colData(ref_sce)$cell_type1 = colData(ref_sce)$Cluster
  rowData(ref_sce)$feature_symbol = rownames(rowData(ref_sce))
  #saveRDS(ref_sce_smillie,"Results-Sam/ref_sce_smillie.RDS")
  # build  the index we will use to map our unlabeled data to. First, we select genes to use, which will be those deemed most informative by scmap after fitting a linear model to the gene expression by gene dropout distribution. Those which are most informative have high expression values and low % dropout rates across cells.
  
  # Create scmap-cluster reference by first selecting the most informative features
  print("selecting the most informative features...")
  ref_sce = scmap::selectFeatures(ref_sce, suppress_plot=FALSE,n_features)
  
  # Create a list of mitochondrial genes from the dataset (genes that begin with "MT")
  print("check id there are mitochondrial genes to filter out...")
  mt_genes = rownames(ref_sce)[grep("^MT-", rownames(ref_sce))]
  # Remove these genes from the features used by scmap
  rowData(ref_sce)$scmap_features[rownames(ref_sce) %in% mt_genes] <- FALSE
  
  # Extract the features and assign them to a new variable, "scmap_feature_genes"
  scmap_feature_genes <- rownames(ref_sce)[which(rowData(ref_sce)$scmap_features)]
  # Note that the number of genes/features is identical to what we just checked

  
  
  if(scmap_cluster){
  # Create reference profiles;
  # Once reference profiles are generated the original data are not needed for scmap-cluster
  print("indexing for scamp cluster...")
  ref_sce <- scmap::indexCluster(ref_sce)
  scmap_cluster_reference <- metadata(ref_sce)$scmap_cluster_index
  }
  
 
  
  if(scmap_cell){
  # Update the previous reference to also contain the scmap-cell reference
  print("indexing for scmap cell...")
  ref_sce <- scmap::indexCell(ref_sce)
  # Extract the scmap index from the reference and store as a variable
  scmap_cell_reference <- metadata(ref_sce)$scmap_cell_index
  # Extract the associated cell IDs from the reference and save as a variable
  scmap_cell_metadata <- colData(ref_sce)
  }
  
  # Now we need to load our unlabeled dataset into R. N
  # Normal preprocessing including QC filtering, normalizing and log-transforming the data must be done prior to annotating.
  # Make SingleCellExperiment from the raw matrix
  print("loading the query object for the annotation...")
  query_sce <- as.SingleCellExperiment(query_seur_obj,assay = query_assay)
  # add feature_symbol column (i.e. the gene symbols)
  rowData(query_sce)$feature_symbol <- rownames(query_sce)
  
  
  if(scmap_cluster){
  # We are now ready to annotate our data with **scmap-cluster**. Let's start with scmap-cluster:
  # Run scmapCluster
  print("Running Scmap cluster...")
  scmap_cluster_res <- scmap::scmapCluster(projection=query_sce, 
                                           index_list=list(scl_ref = scmap_cluster_reference), 
                                           threshold=0.1)
  # Store this annotation information within the query object
  colData(query_sce)$scmap_cluster <- scmap_cluster_res$combined_labs
  query_sce$scmap_cluster = gsub("ambiguous", "unknown", query_sce$scmap_cluster)
  query_sce$scmap_cluster = gsub("unassigned", "unknown", query_sce$scmap_cluster)
  query_seur_obj$scmap_cluster = query_sce$scmap_cluster
  }
  
  
  if(scmap_cell){
  #Alternatively we could use scmap-cell, to find the 10 nearest neighbours to each cell (i.e. the 10 most similar cells to each query cell), 
  # then pick the annotation that is most common among the neighbours:
  # Determine the 10 nearest neighbours from the reference dataset for each
  # cell in the query dataset using scmapCell
  nearest_neighbours <- scmap::scmapCell(projection=query_sce, 
                                         index_list = list(sce_ref = scmap_cell_reference), 
                                         w=10)
  # Get metadata (cell type IDs) for the neighbours of each cell in the query dataset
  mode_label <- function(neighbours, metadata=scmap_cell_metadata$cell_type1) {
    freq <- table(metadata[neighbours])
    label <- names(freq)[which(freq == max(freq))]
    if (length(label) > 1) {return("ambiguous")}
    return(label)
  }
  # Apply these labels to the query cells
  scmap_cell_labs <- apply(nearest_neighbours$sce_ref$cells, 2, mode_label)
  
  # Add the labels to the query object
  colData(query_sce)$scmap_cell <- scmap_cell_labs
  query_sce$scmap_cell = gsub("ambiguous", "unknown", query_sce$scmap_cell)
  query_sce$scmap_cell = gsub("unassigned", "unknown", query_sce$scmap_cell)
  query_seur_obj$scmap_cell = query_sce$scmap_cell
  }
  
  return(query_seur_obj)
  
}



#### assigning input parameters properly based on the version of SingleR which needs to be run (you can run one of the SingleR versions by setting others to False):
# 1 - Default version of SingleR : Set to False (default=False) if you want to skip this run
# 2 - SingleR with input list of celltype-specific markers per celltype: Set to False (add_markers=False) if you want to skip this version, otherwise make sure to assign "marker_list" based on the example file 
# 3 - SingleDE with input list of DE markers: Setto False (add_DEGs= False) if you want to skip this version, otherwise make sure to assign "DEGs_list" based on the example file
# please make sure at least one of the SingleR versions is properly set to be run...

# Cluster: is a name of the column in reference Seurat object which includes the reference cell labels


Run_SingleR = function (ref_seur_obj, query_seur_obj, ref_assay ="RNA", 
                        query_assay = "RNA", default=T ,add_markers = T,
                        marker_list, add_DEGs = T, DEGs_list, Cluster,
                        cellTypes_2_exclude = c("Glia")){
  
  print("loading the reference object for the annotation...")
  ref_sce = as.SingleCellExperiment(ref_seur_obj,ref_assay)
  #colnames(colData(ref_sce))
  colData(ref_sce)$Cluster = colData(ref_sce)[[Cluster]]
  
  # Assign cell-type labels in a column named "cell_type1"
  colData(ref_sce)$cell_type1 = colData(ref_sce)$Cluster
  rowData(ref_sce)$feature_symbol = rownames(rowData(ref_sce))
  
  
  print("loading the query object for the annotation...")
  query_sce <- as.SingleCellExperiment(query_seur_obj,assay = query_assay)
  # add feature_symbol column (i.e. the gene symbols)
  rowData(query_sce)$feature_symbol <- rownames(query_sce)
  
  
  
  if(default){
    
    print("Running SingleR wihtout including marker genes...")
    # Run SingleR on the query data and the reference to acquire
    # cell-type predictions for the cells in the query dataset
    predictions <- SingleR::SingleR(test=query_sce, ref=ref_sce, labels=ref_sce$cell_type1)
    # You'll notice that some of the cells didn't get assigned a cell identity
    # Change NAs to "ambiguous"
    predictions$pruned.labels[which(is.na(predictions$pruned.labels))] <- "unknown"
    # Add singleR labels to query_sce
    colData(query_sce)$singleR <- predictions$pruned.labels
    query_seur_obj$singleR = query_sce$singleR
  }
  
  if(add_markers){
    print("Running SingleR with celltype-specific markers...")
    print("SingleR training with markers...")
    singleR_trained = trainSingleR(
      ref=ref_sce,
      labels=ref_sce$cell_type1,
      genes = marker_list
    )
    
    print("Query data classification using SingleR...")
    predictions = classifySingleR(query_sce, singleR_trained) 
    # You'll notice that some of the cells didn't get assigned a cell identity
    # We can count the number here:
    sum(is.na(predictions$pruned.labels))
    # Change NAs to "ambiguous"
    predictions$pruned.labels[which(is.na(predictions$pruned.labels))] <- "unknown"
    # Add singleR labels to query_sce
    colData(query_sce)$singleRm <- predictions$pruned.labels
    query_seur_obj$singleRm = query_sce$singleRm
  }
  
  if(add_DEGs){
    print("Running SingleR with pairwise DE markers...")
    SingleR.pred = SingleR_annotation(query_seur_obj, ref_seur_obj,
                                      cellTypes_2_exclude, 
                                      DEGs_list)
    colData(query_sce)$singleRDE <- SingleR.pred$labels
    query_seur_obj$singleRDE = query_sce$singleRDE
    
  }
  
  return(query_seur_obj)

}

#the labels should be under the Cluster in the reference object (ref_seur_obj)
SingleR_annotation <- function(query_seur_obj, # New dataset to annotate
                               ref_seur_obj, # Reference set
                               cellTypes_2_exclude = c("Glia"),
                               DEGs_list, # Selected DE Markers to refine results from cell-cell pairwise T-test
                               #Note that DEGs_list here is in the format of DE analysis/ MAST output
                               N.top.markers = 50,
                               quantile=0.95){
  
  # transform query_seur_obj to SingleCellExperiment
  so.sce <- as.SingleCellExperiment(query_seur_obj,assay = "RNA")
  
  # Training set - in this case Smillie et al
  train.sce <- as.SingleCellExperiment(ref_seur_obj[,!(ref_seur_obj$Cluster %in% cellTypes_2_exclude)],assay = "RNA" )
  
  # For the SingleR run, retain only genes that are shared between training (reference) and test (Vedo)
  common <- intersect(rownames(train.sce), rownames(so.sce))
  train.sce <- train.sce[common,]
  so.sce <- so.sce[common,]
  
  
  # List 1 - Get top 50 markers per Cluster in the ref_seur_obj
  #perform pairwise TTest - can be changed if needed to another test
  out <- pairwiseTTests(logcounts(train.sce), train.sce$Cluster, direction="up")
  markers.list <- getTopMarkers(out$statistics, out$pairs, n=N.top.markers)
  
  # List 2- Create a DE.gene list based on DE results from Smillie et al.or in-house DE analysis MAST etc 
  de.genes <- list()
  DEGs_list <- DEGs_list[DEGs_list$ident %in% unique(Train.seur$Cluster),]
  for(i in unique(DEGs_list$ident)){
    de.genes[[i]] <- DEGs_list[DEGs_list$ident ==i,]$gene
  }
  

  # Combine List 1 & List 2
  
  for(i in c(names(de.genes))){
    cell.type <- i
    cell.type.markers <-  de.genes[[cell.type]]
    markers.list[[cell.type]] <- lapply(seq_along(names(markers.list[[cell.type]])), function(idx) {cell.type.markers[(cell.type.markers %in% markers.list[[cell.type]][[idx]]) & cell.type.markers %in% common]})
    names(markers.list[[cell.type]]) <- names(de.genes)
  }
  
  #trainSingleR for the Reference Smilie Fib dataset
  trained <- trainSingleR(train.sce, labels=train.sce$Cluster, genes=markers.list)
  #PredictSingleR for the Test Vedolizumab Fib dataset
  SingleR.pred = classifySingleR(so.sce, trained,quantile=quantile) 
  
  
  return(SingleR.pred)
  
}



get_consensus_label <- function(labels){
  labels <- labels[labels != "unknown"]
  if (length(labels) == 0) {return("unknown")} # No annotation from the methods
  freq <- table(labels)
  label <- names(freq)[which(freq == max(freq))]
  if (length(label) > 1) {return("unknown")} # if there is no winner among the annotations
  return(label)
}
