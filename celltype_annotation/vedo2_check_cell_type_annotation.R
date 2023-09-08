######################
# libraries          #
######################

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)


####################
# Functions        #
####################

create_confusion_matrix <- function(assignment_table, truth_column, prediction_column, truth_column_label=NULL, prediction_column_label=NULL){
  # init the table
  confusion_table <- NULL
  # check each truth
  for(truth in unique(assignment_table[[truth_column]])){
    # get these truths
    truth_rows <- assignment_table[assignment_table[[truth_column]] == truth, ]
    # check now many have this truth
    this_truth_number <- nrow(truth_rows)
    # check what was predicted for these truths
    #for(prediction in unique(truth_rows[[prediction_column]])){
    for(prediction in unique(assignment_table[[prediction_column]])){
      # check the number of this prediction
      this_prediction_number <- nrow(truth_rows[truth_rows[[prediction_column]] == prediction, ])
      # init variable
      fraction <- NULL
      # we can only calculate a fraction if the result is not zero
      if(this_prediction_number > 0){
        # calculate the fraction
        fraction <- this_prediction_number / this_truth_number
      }
      # otherwise we just set it to zero
      else{
        fraction <- 0
      }
      # turn into row
      this_row <- data.frame(truth=c(truth), prediction=c(prediction), freq=c(fraction), stringsAsFactors = F)
      # add this entry to the dataframe
      if(is.null(confusion_table)){
        confusion_table <- this_row
      }
      else{
        confusion_table <- rbind(confusion_table, this_row)
      }
    }
  }
  # round the frequency off to a sensible cutoff
  confusion_table$freq <- round(confusion_table$freq, digits=2)
  # turn into plot
  p <- ggplot(data=confusion_table, aes(x=truth, y=prediction, fill=freq)) + geom_tile() + scale_fill_gradient(low='red', high='blue') + geom_text(aes(label=freq))
  # some options
  if(!is.null(truth_column_label)){
    p <- p + xlab(truth_column_label)
  }
  if(!is.null(prediction_column_label)){
    p <- p + ylab(prediction_column_label)
  }
  return(p)
}

####################
# Main Code        #
####################

# storage spots
storage_read_part <- 'tmp01'
storage_write_part <- 'tmp01'

# object locations
object_loc <- paste('/groups/umcg-weersma/', storage_read_part, '/projects/vedo2_predict/ongoing/seurat_preprocess_samples/objects/', sep = '')
vedo2_object_loc <- paste(object_loc, 'vedo2_batch1_biop_integrated_noribo_nonorm_withmito_ref.rds', sep = '')


metadata <- read.csv('~/Desktop/vedo2_takeda_UMCG_ct_classification.csv')
