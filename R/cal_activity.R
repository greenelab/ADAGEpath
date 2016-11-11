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


#' Signature extraction
#'
#' Extracts all gene signatures from an ADAGE model.
#'
#' @param model the ADAGE model to be used for extracting gene signatures
#' (default: the 300-node eADAGE model preloaded in the pacakge).
#' @param HW_cutoff number of standard deviations from mean in a node's weight
#' distribution to be considered as high-weight (default to 2.5). Only
#' high-weight genes are included in gene signatures.
#' @param use_symbol logical, whether the returned signatures use gene symbol
#' as gene identifiers.
#' @return a named list with each element being a gene signature
extract_signatures <- function(model = eADAGEmodel, HW_cutoff = 2.5,
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


#' One signature extraction
#'
#' Extracts a single gene signature from an ADAGE model.
#'
#' @param node_weight a vector storing weight values of each gene to a node
#' @param geneID gene identifiers correspond to the weight vector
#' @param side character, "pos" or "neg"
#' @param HW_cutoff number of standard deviations from mean in a node's weight
#' distribution to be considered as high-weight (default to 2.5).
#' @return a character vector storing genes in the signature defined by the
#' weight vector and the side.
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


enrich_test <- function(set1, set2, set_all){

  high_in <- length(intersect(set1, set2))
  low_in <- length(set1) - high_in
  high_out <- length(set2) - high_in
  low_out <- length(set_all) - high_in - high_out - low_in
  contingency_table <- matrix(c(high_in, low_in, high_out, low_out), nrow = 2)
  pvalue <- fisher.test(contingency_table, alternative='greater')$p.value

  return(pvalue)
}

calculate_signature_similarity <- function(signature_list){
  comb_sig <- combn(signature_list, 2)
  all_HWGs <- unique(unlist(signature_list))
  # TODO: use parellel version of mapply
  pvalue_list <- mapply(enrich_test, comb_sig[1, ], comb_sig[2, ],
                        MoreArgs = list(set_all = all_HWGs))
  qvalue_list <- p.adjust(pvalue_list, method = "fdr")

  return(qvalue_list)

}


build_signature_similarity_matrix <- function(signature_list, overlap_qvalues){
  overlapQ_matrix = matrix(, nrow = length(signature_list),
                           ncol = length(signature_list))
  overlapQ_matrix[lower.tri(overlapQ_matrix, diag = FALSE)] <- overlap_qvalues
  overlapQ_matrix[upper.tri(overlapQ_matrix)] <-
    t(overlapQ_matrix)[upper.tri(overlapQ_matrix)]
  diag(overlapQ_matrix) <- min(overlapQ_matrix, na.rm = TRUE)
  rownames(overlapQ_matrix) <- names(signature_list)
  colnames(overlapQ_matrix) <- names(signature_list)

  return(overlapQ_matrix)
}


plot_signature_similarity <- function(overlapQ_matrix, signatures = NULL){
  if (!is.null(signatures)){
    if (signatures %in% rownames(overlapQ_matrix)){
      overlapQ_matrix <- overlapQ_matrix[signatures, signatures]
    } else {
      stop("Given signatures are not in the overlap matrix!")
    }

  }

  overlapQ_matrix <- -log10(overlapQ_matrix)
  corrplot::corrplot(overlapQ_matrix, is.corr = FALSE, method = 'square',
                     order = 'hclust', diag = FALSE)
}
