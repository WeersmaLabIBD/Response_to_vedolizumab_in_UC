#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Paola Pibiri, Roy Oelen
# Name: vedo2_nichenet_lr_circles.R
# Function: create nichenet circleplots for vedo2
############################################################################################################################


####################
# libraries        #
####################

library(nichenetr)
library(Seurat)
library(tidyverse)
library(circlize)
library(RColorBrewer)


####################
# functions        #
####################


get_link_info_receiver <- function(ligand_target_matrix, lr_network, weighted_networks,
                                      seuratObj, 
                                      receiver, 
                                      senders=c("Post-capillary Venules","IgG", "M cells", "NKs", "CD8+ IL17+", "DCs", "ILCs"), 
                                      condition_colname='condition', 
                                      condition_oi='no_T0',
                                      condition_reference='yes_T0',
                                      top_n_ligands=30
                                      ) {
  # predict active target genes and construct an active ligand-receptor network
  nichenet_output = nichenet_seuratobj_aggregate(
    seurat_obj = seuratObj, 
    receiver = receiver, 
    condition_colname = condition_colname, condition_oi = condition_oi, condition_reference = condition_reference, 
    sender = sender_celltypes, 
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks,
    top_n_ligands = top_n_ligands, 
    assay_oi = 'SCT')
  
  # get the average expression of the ligands per cell type
  avg_expression_ligands = AverageExpression(seuratObj, features = nichenet_output$top_ligands)
  
  # look for which sender cell types show an expression that is higher than the average + SD
  sender_ligand_assignment = avg_expression_ligands$SCT %>% apply(1, function(ligand_expression){
    ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
  }) %>% t()
  sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
  
  # get all ligands
  all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
  # as well as the unique ones
  unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
  # and the general ones
  general_ligands = nichenet_output$top_ligands %>% setdiff(unique_ligands)
  
  # put the results in a list
  result_list <- list(
    'nichenet_output' = nichenet_output,
    'avg_expression_ligands' = avg_expression_ligands,
    'sender_ligand_assignment' = sender_ligand_assignment,
    'all_assigned_ligands' = all_assigned_ligands,
    'unique_ligands' = unique_ligands,
    'general_ligands' = general_ligands
  )
  return(result_list)
}


get_link_info_per_receiver <- function(ligand_target_matrix, lr_network, weighted_networks,
                                       seuratObj, 
                                       receivers=c("Inflammatory Monocytes", "Inflammatory Fibroblasts", "Tregs"), 
                                       senders=c("Post-capillary Venules","IgG", "M cells", "NKs", "CD8+ IL17+", "DCs", "ILCs"), 
                                       condition_colname='condition', 
                                       condition_oi='no_T0',
                                       condition_reference='yes_T0',
                                       top_n_ligands=30
) {
  # save the result per receiver
  result_per_receiver <- list()
  # check each receiver
  for (receiver in receivers) {
    tryCatch(
      {
        result_receiver <- get_link_info_receiver(ligand_target_matrix, lr_network, weighted_networks,
                                                 seuratObj=seuratObj, 
                                                 receiver=receiver, 
                                                 senders=senders, 
                                                 condition_colname=condition_colname, 
                                                 condition_oi=condition_oi,
                                                 condition_reference=condition_reference,
                                                 top_n_ligands=top_n_ligands
        )
        # add the result to the list
        result_per_receiver[[receiver]] <- result_receiver
      }, error=function(cond){
        print(paste0('problem with receiver: ', receiver))
        message(cond)
      }
    )
        
    }
    return(result_per_receiver)
}

add_sender_info <- function(ligand_target_df, sender_ligand_assignment) {
  # make a list of dataframes with the relevant cell type and ligands
  ligands_per_ct_df_list <- list()
  # check each cell type
  for (cell_type in names(sender_ligand_assignment)) {
    # extract the ligands
    ligands <- sender_ligand_assignment[[cell_type]]
    # get the ones which are true
    ligands <- ligands[ligands == T]
    # put this in a dataframe
    ligands_ct <- data.frame('ligand_type' = cell_type, 'ligand' = names(ligands))
    # add to the list
    ligands_per_ct_df_list[[cell_type]] <- ligands_ct
  }
  # add all together
  ligands_per_ct_df <- do.call('rbind',ligands_per_ct_df_list)
  ligand_target_with_type <- ligand_target_df %>% inner_join(ligands_per_ct_df)
  return(ligand_target_with_type)
}

add_sender_info_per_receiver <- function(link_info_per_receiver) {
  # check each receiver
  for (receiver in names(link_info_per_receiver)) {
    # extract for this receiver
    info_receiver <- link_info_per_receiver[[receiver]]
    # add the target type, which is the receiver
    print(nrow(info_receiver[['nichenet_output']]$ligand_target_df))
    active_ligand_target_links_df_receiver = info_receiver[['nichenet_output']]$ligand_target_df %>% mutate(target_type = receiver)
    # and add the ligand type, which was the sending cell
    active_ligand_target_links_df_receiver <- add_sender_info(active_ligand_target_links_df_receiver, info_receiver[['sender_ligand_assignment']])
    # we can do the same with the receptors of course
    active_ligand_receptor_links_df_receiver = info_receiver[['nichenet_output']]$ligand_receptor_df %>% mutate(target_type = receiver)
    active_ligand_receptor_links_df_receiver <- add_sender_info(active_ligand_receptor_links_df_receiver, info_receiver[['sender_ligand_assignment']])
    # add as new values
    info_receiver[['ligand_target_df_wreceiver']] <- active_ligand_target_links_df_receiver
    info_receiver[['ligand_receptor_df_wreceiver']] <- active_ligand_receptor_links_df_receiver
    # add back into the list
    link_info_per_receiver[[receiver]] <- info_receiver
  }
  return(link_info_per_receiver)
}

merging_info <- function(link_info_per_receiver) {
  list_df <- list()
  for (celltype in names(link_info_per_receiver)) {
    link_info_receiver <- link_info_per_receiver[[celltype]][['ligand_receptor_df_wreceiver']]
    list_df[[celltype]] <- link_info_receiver
  }
  merged <- do.call('rbind', list_df)
  return(merged)
}


info_to_plot_table <- function(info_table){
  # make into the format specified by nichenet
  plot_table <- data.frame(
    'sender' = info_table$ligand_type,
    'receiver' = info_table$target_type,
    'ligand' = info_table$ligand,
    'receptor' = info_table$receptor,
    'ligand_receptor' = paste(info_table$ligand, info_table$receptor, sep = "--"),
    'niche' = paste(info_table$ligand_type, "niche", sep = "_"),
    'prioritization_score' = info_table$weight,
    'top_niche' = paste(info_table$ligand_type, "niche", sep = "_")
  )
  return(plot_table)
}

####################
#  main code       #
####################

# location of object
vedo2_object_loc <- '/groups/umcg-weersma/tmp01/projects/vedopredict2/biopsies/ongoing/seurat_preprocess_samples/objects/vedo2_remerged_resct_20230215.rds'


vedo2 <- readRDS(vedo2_object_loc)

# create the new celltypes
vedo2@meta.data$celltype_final.v1 <- vedo2@meta.data$celltype_final
vedo2@meta.data[vedo2@meta.data$celltype_final.v1 %in% c('DC1', 'DC2'), 'celltype_final.v1']  <- 'DCs'

# add a column denoting both the response and the timepoint
vedo2@meta.data$condition <- paste0(vedo2@meta.data$PGA._resp, '_', vedo2@meta.data$timepoint)

# use the cell types as the identities
Idents(object = vedo2) <- 'celltype_final.v1'

# these are our sender cell types
sender_celltypes = c("Post-capillary Venules","IgG", "M cells", "NKs", "CD8+ IL17+", "DCs", "ILCs")
# and these are the receiver cell types
receiver = c("Inflammatory Monocytes", "Inflammatory Fibroblasts", "Tregs")

#subset the object for only my relevant celltypes, should be the join of the senders and the receivers
relevant_celltypes <- c('Inflammatory Monocytes', 'NKs', 'Tregs', 'Inflammatory Fibroblasts', 'Post-capillary Venules', 'DCs', 'IgG', 'M cells', 'CD8+ IL17+', 'ILCs')
vedo2 <- vedo2[ ,vedo2@meta.data$celltype_final.v1 %in% relevant_celltypes]


# read the annotations of Nichenet
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))




# # test with tregs as receiver
# info_tregs <- get_link_info_receiver(ligand_target_matrix = ligand_target_matrix,
#                        lr_network = lr_network,
#                        weighted_networks = weighted_networks,
#                        seuratObj = vedo2,
#                        receiver = 'Tregs',
#                        top_n_ligands = 10)
# 
# # add the target type, which is the receiver
# active_ligand_target_links_df_tregs = info_tregs[['nichenet_output']]$ligand_target_df %>% mutate(target_type = "Tregs")
# # and add the ligand type, which was the sending cell
# active_ligand_target_links_df_tregs <- add_sender_info(active_ligand_target_links_df_tregs, info_tregs[['sender_ligand_assignment']])
# # we can do the same with the receptors of course
# active_ligand_receptor_links_df_tregs = info_tregs[['nichenet_output']]$ligand_receptor_df %>% mutate(receptor_type = "Tregs")
# active_ligand_receptor_links_df_tregs <- add_sender_info(active_ligand_receptor_links_df_tregs, info_tregs[['sender_ligand_assignment']])
# 
# 
# #repeat for inflammatory monocytes 
# 
# # test with tregs as receiver
# info_mono <- get_link_info_receiver(ligand_target_matrix = ligand_target_matrix,
#                                      lr_network = lr_network,
#                                      weighted_networks = weighted_networks,
#                                      seuratObj = vedo2,
#                                      receiver = 'Inflammatory Monocytes',
#                                      top_n_ligands = 10)
# # add the target type, which is the receiver
# active_ligand_target_links_df_mono = info_mono[['nichenet_output']]$ligand_target_df %>% mutate(target_type = "Inflammatory Monocytes")
# # and add the ligand type, which was the sending cell
# active_ligand_target_links_df_mono <- add_sender_info(active_ligand_target_links_df_mono, info_mono[['sender_ligand_assignment']])
# # we can do the same with the receptors of course
# active_ligand_receptor_links_df_mono = info_mono[['nichenet_output']]$ligand_receptor_df %>% mutate(receptor_type = "Inflammatory Monocytes")
# active_ligand_receptor_links_df_mono <- add_sender_info(active_ligand_receptor_links_df_mono, info_mono[['sender_ligand_assignment']])
# 
# 
# #same for inflammatory fibroblasts
# 
# # test with tregs as receiver
# info_fibro <- get_link_info_receiver(ligand_target_matrix = ligand_target_matrix,
#                                     lr_network = lr_network,
#                                     weighted_networks = weighted_networks,
#                                     seuratObj = vedo2,
#                                     receiver = 'Inflammatory Fibroblasts',
#                                     top_n_ligands = 10)
# # add the target type, which is the receiver
# active_ligand_target_links_df_fibro = info_fibro[['nichenet_output']]$ligand_target_df %>% mutate(target_type = "Inflammatory Fibroblasts")
# # and add the ligand type, which was the sending cell
# active_ligand_target_links_df_fibro <- add_sender_info(active_ligand_target_links_df_fibro, info_fibro[['sender_ligand_assignment']])
# # we can do the same with the receptors of course
# active_ligand_receptor_links_df_fibro = info_fibro[['nichenet_output']]$ligand_receptor_df %>% mutate(receptor_type = "Inflammatory Fibroblasts")
# active_ligand_receptor_links_df_fibro <- add_sender_info(active_ligand_receptor_links_df_fibro, info_fibro[['sender_ligand_assignment']])
# 
# #merge the 3 tables together
# merged <- do.call('rbind', list(active_ligand_receptor_links_df_tregs, active_ligand_receptor_links_df_mono, active_ligand_receptor_links_df_fibro))
# 
# # make into the format specified by nichenet
# t0_NR_vs_R <- data.frame(
#   'sender' = merged$ligand_type,
#   'receiver' = merged$receptor_type,
#   'ligand' = merged$ligand,
#   'receptor' = merged$receptor,
#   'ligand_receptor' = paste(merged$ligand, merged$receptor, sep = "--"),
#   'niche' = paste(merged$ligand_type, "niche", sep = "_"),
#   'prioritization_score' = merged$weight,
#   'top_niche' = paste(merged$ligand_type, "niche", sep = "_")
# )
# 
# t0_NR_vs_R <- t0_NR_vs_R[!(t0_NR_vs_R$sender %in% receiver), ]
# 
# # create the colour palette
# 
# colors_celltypes <- (
# brewer.pal(
# n = length(unique(c(t0_NR_vs_R$sender, t0_NR_vs_R$receiver))), name = 'Paired')
# )
# names(colors_celltypes) <- unique(c(t0_NR_vs_R$sender, t0_NR_vs_R$receiver))
# 
# # sort by weight
# t0_NR_vs_R[order(t0_NR_vs_R$prioritization_score, decreasing = T), ]
# # now remove duplicates, keeping the first one (which is the strongest due to sorting)
# df <- df[!duplicated(df$ligand_receptor), ]
# 
# # make the plot
# par(mfrow = c(1, 2))
# make_circos_lr(df, colors_celltypes, colors_celltypes, transparency = 0.3, scale = 0.1)
# par(mfrow = c(1, 1))
# 



####################################################### automate everything

# NR vs R at baseline, core cells as receivers 
sender_celltypes = c("Post-capillary Venules","IgG", "M cells", "NKs", "CD8+ IL17+", "DCs", "ILCs")
# and these are the receiver cell types
receiver = c("Inflammatory Monocytes", "Inflammatory Fibroblasts", "Tregs")

# get some colours
colors_celltypes <- (
  brewer.pal(
    n = length(unique(c(sender_celltypes, receiver))), name = 'Paired'
  )
)
names(colors_celltypes) <- unique(c(sender_celltypes, receiver))


# get prioritised ligands and receptors per receiver
t0_NR_vs_R_per_receiver <- get_link_info_per_receiver(ligand_target_matrix, lr_network, weighted_networks, vedo2, 
                                                      receivers = receiver, senders = senders, 
                                                      top_n_ligands = 10)
# add extra info
t0_NR_vs_R_per_receiver_infoadded <- add_sender_info_per_receiver(t0_NR_vs_R_per_receiver)
# merging
t0_NR_vs_R_merged <- merging_info(t0_NR_vs_R_per_receiver_infoadded)
# make plottable
t0_NR_vs_R_plottable <- info_to_plot_table(t0_NR_vs_R_merged)
# remove the cell types that we only want as receiver
t0_NR_vs_R_plottable <- t0_NR_vs_R_plottable[!(t0_NR_vs_R_plottable$sender %in% receiver), ]
# sort by weight
t0_NR_vs_R_plottable[order(t0_NR_vs_R_plottable$prioritization_score, decreasing = T), ]
# now remove duplicates, keeping the first one (which is the strongest due to sorting)
t0_NR_vs_R_plottable <- t0_NR_vs_R_plottable[!duplicated(t0_NR_vs_R_plottable$ligand_receptor), ]
# make the plot
par(mfrow = c(1, 2))
make_circos_lr(t0_NR_vs_R_plottable, colors_celltypes, colors_celltypes, transparency = 0.3, scale = 0.1)
par(mfrow = c(1, 1))


#NR vs R at baseline, core cells as senders
receiver_rev = c("Post-capillary Venules","IgG", "M cells", "NKs", "CD8+ IL17+", "DCs", "ILCs")
# and these are the receiver cell types
sender_celltypes_rev = c("Inflammatory Monocytes", "Inflammatory Fibroblasts", "Tregs")
# do all the steps again
t0_NR_vs_R_per_receiver_rev <- get_link_info_per_receiver(ligand_target_matrix, lr_network, weighted_networks, vedo2, 
                                                      receivers = receiver_rev, senders = senders_rev, 
                                                      top_n_ligands = 5)
t0_NR_vs_R_per_receiver_infoadded_rev <- add_sender_info_per_receiver(t0_NR_vs_R_per_receiver_rev)
t0_NR_vs_R_merged_rev <- merging_info(t0_NR_vs_R_per_receiver_infoadded_rev)
t0_NR_vs_R_plottable_rev <- info_to_plot_table(t0_NR_vs_R_merged_rev)
t0_NR_vs_R_plottable_rev <- t0_NR_vs_R_plottable_rev[!(t0_NR_vs_R_plottable_rev$sender %in% receiver_rev), ]
t0_NR_vs_R_plottable_rev[order(t0_NR_vs_R_plottable_rev$prioritization_score, decreasing = T), ]
t0_NR_vs_R_plottable_rev <- t0_NR_vs_R_plottable_rev[!duplicated(t0_NR_vs_R_plottable_rev$ligand_receptor), ]
# make the plot
par(mfrow = c(1, 2))
make_circos_lr(t0_NR_vs_R_plottable_rev, colors_celltypes, colors_celltypes, transparency = 0.3, scale = 0.1)
par(mfrow = c(1, 1))



# responders T4 vs T0
R_T4_vs_T0_per_receiver <- get_link_info_per_receiver(ligand_target_matrix, lr_network, weighted_networks, vedo2, 
                                                      receivers = receiver, senders = senders, 
                                                      top_n_ligands = 10,
                                                      condition_oi='yes_T4',
                                                      condition_reference='yes_T0')
R_T4_vs_T0_per_receiver_infoadded <- add_sender_info_per_receiver(R_T4_vs_T0_per_receiver)
R_T4_vs_T0_merged <- merging_info(R_T4_vs_T0_per_receiver_infoadded)
R_T4_vs_T0_plottable <- info_to_plot_table(R_T4_vs_T0_merged)
R_T4_vs_T0_plottable <- R_T4_vs_T0_plottable[!(R_T4_vs_T0_plottable$sender %in% receiver), ]
R_T4_vs_T0_plottable[order(R_T4_vs_T0_plottable$prioritization_score, decreasing = T), ]
R_T4_vs_T0_plottable <- R_T4_vs_T0_plottable[!duplicated(R_T4_vs_T0_plottable$ligand_receptor), ]
# make the plot
par(mfrow = c(1, 2))
make_circos_lr(R_T4_vs_T0_plottable, colors_celltypes, colors_celltypes, transparency = 0.3, scale = 0.1)
par(mfrow = c(1, 1))


# non-reponders T4 vs T0
NR_T4_vs_T0_per_receiver <- get_link_info_per_receiver(ligand_target_matrix, lr_network, weighted_networks, vedo2, 
                                                      receivers = receiver, senders = senders, 
                                                      top_n_ligands = 5,
                                                      condition_oi='no_T4',
                                                      condition_reference='no_T0')
NR_T4_vs_T0_per_receiver_infoadded <- add_sender_info_per_receiver(NR_T4_vs_T0_per_receiver)
NR_T4_vs_T0_merged <- merging_info(NR_T4_vs_T0_per_receiver_infoadded)
NR_T4_vs_T0_plottable <- info_to_plot_table(NR_T4_vs_T0_merged)
NR_T4_vs_T0_plottable <- NR_T4_vs_T0_plottable[!(NR_T4_vs_T0_plottable$sender %in% receiver), ]
NR_T4_vs_T0_plottable[order(NR_T4_vs_T0_plottable$prioritization_score, decreasing = T), ]
NR_T4_vs_T0_plottable <- NR_T4_vs_T0_plottable[!duplicated(NR_T4_vs_T0_plottable$ligand_receptor), ]
# make the plot
par(mfrow = c(1, 2))
make_circos_lr(NR_T4_vs_T0_plottable, colors_celltypes, colors_celltypes, transparency = 0.3, scale = 0.1)
par(mfrow = c(1, 1))