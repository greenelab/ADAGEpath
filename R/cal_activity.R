#' Signature activity calculation
#'
#' Calculates activities for each signature in an ADAGE model
#' specified by the weight matrix for each sample in the input data
#'
#' @param input_data a data frame with gene IDs in the first column and
#' expression values from the second column.
#' @param model the ADAGE model to be used for calculating signature activity
#' (default: the 300-node eADAGE model preloaded in the pacakge).
#' @param HW_cutoff number of standard deviations from mean in a node's weight
#' distribution to be considered as high-weight (default to 2.5).
#' @return a named matrix storing activities per signature per sample
#' @export
calculate_activity <- function(input_data, model = eADAGEmodel,
                               HW_cutoff = 2.5) {

  if (!check_input(input_data)){
    stop("The input data should be a data frame with first column storing
         geneIDs in character and the rest columns storing expression values
         for each sample in numeric.")
  }

  if(!all(input_data[, 1] == model[, 1])){
    stop("The gene IDs stored in the first column of the input data and the
         first column of the model should be the same.")
  }

  weight_matrix <- as.matrix(model[, -1])
  model_size <- ncol(weight_matrix)
  value_only <- as.matrix(input_data[, -1])

  HWactivity_perGene_pos <- sapply(1:model_size, function(x)
    one_signature_activity(weight_matrix = weight_matrix,
                           express_matrix = value_only,
                           node = x, side = "pos", HW_cutoff = HW_cutoff))
  HWactivity_perGene_neg <- sapply(1:model_size, function(x)
    one_signature_activity(weight_matrix = weight_matrix,
                           express_matrix = value_only,
                           node = x, side = "neg", HW_cutoff = HW_cutoff))

  # combine positive and negative sides
  HWactivity_perGene <- dplyr::bind_cols(data.frame(HWactivity_perGene_pos),
                                         data.frame(HWactivity_perGene_neg))
  rownames(HWactivity_perGene) <- colnames(value_only)
  colnames(HWactivity_perGene) <- c(paste0("Node", seq(1, model_size), "pos"),
                                    paste0("Node", seq(1, model_size), "neg"))

  # omit positive and negative signs
  HWactivity_perGene <- abs(HWactivity_perGene)

  return(t(HWactivity_perGene))

}


#' One signature activity
#'
#' Calculates activities for a specific signature in an ADAGE model
#'
#' @param weight_matrix a data matrix storing the weight matrix in an ADAGE model.
#' @param express_matrix a data matrix storing gene expression values in a dataset
#' @param node a int ranging from 1 to number of columns in weight_matrix
#' @param side character, "pos" or "neg"
#' @param HW_cutoff number of standard deviations from mean in a node's weight
#' distribution to be considered as high-weight (default to 2.5).
#' @return a vector storing activities of one signature across samples
one_signature_activity <- function(weight_matrix, express_matrix, node, side,
                                   HW_cutoff = 2.5){

  if (node > ncol(weight_matrix)){
    stop("Node too large, no such node in the provided weight matrix!")
  }

  if (side == "pos") {
    pos_cutoff <- mean(weight_matrix[, node]) + HW_cutoff *
      sd(weight_matrix[, node])
    HWG_index <- weight_matrix[, node] >= pos_cutoff
  } else if (side == "neg"){
    neg_cutoff <- mean(weight_matrix[, node]) - HW_cutoff *
      sd(weight_matrix[, node])
    HWG_index <- weight_matrix[, node] <= neg_cutoff
  }

  node_activity <- t(weight_matrix[HWG_index, node]) %*%
    express_matrix[HWG_index, ]
  if (any(HWG_index) > 0) {
    node_activity <- node_activity / sum(HWG_index)
  }

  return(node_activity)
}


plot_activity_heatmap <- function(){

}
