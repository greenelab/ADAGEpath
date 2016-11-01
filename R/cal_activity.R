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


get_signatures <- function(model = eADAGEmodel, HW_cutoff = 2.5,
                           use_symbol = FALSE){

  geneID <- as.data.frame(model)[, 1]
  if (use_symbol) {
    geneID <- sapply(geneID, function(x) to_symbol(x))
  }

  weight_matrix <- as.matrix(model[, -1])
  model_size <- ncol(weight_matrix)

  pos_signatures <- lapply(1:model_size, function(x)
    one_signature(weight_matrix[, x], geneID, "pos", HW_cutoff))
  neg_signatures <- lapply(1:model_size, function(x)
    one_signature(weight_matrix[, x], geneID, "neg", HW_cutoff))

  signature_list <- c(pos_signatures, neg_signatures)
  names(signature_list) <- c(paste0("Node", seq(1, model_size), "pos"),
                             paste0("Node", seq(1, model_size), "neg"))

  return(signature_list)
}


one_signature <- function(node_weight, geneID, side, HW_cutoff = 2.5){

  if (side == "pos") {

    pos_cutoff <- mean(node_weight) + HW_cutoff * sd(node_weight)
    # order positive HW genes from high to low weight
    HW_order <- order(node_weight, decreasing = TRUE)
    ordered_node_weight <- node_weight[HW_order]
    ordered_geneID <- geneID[HW_order]
    HWG <- ordered_geneID[ordered_node_weight >= pos_cutoff]

  } else if (side == "neg") {

    neg_cutoff <- mean(node_weight) - HW_cutoff * sd(node_weight)
    # order negative HW genes from low to high weight
    HW_order <- order(node_weight, decreasing = FALSE)
    ordered_node_weight <- node_weight[HW_order]
    ordered_geneID <- geneID[HW_order]
    HWG <- ordered_geneID[ordered_node_weight <= neg_cutoff]

  } else {
    stop("side can only be pos or neg.")
  }

  return(HWG)
}


