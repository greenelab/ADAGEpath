
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
cal_activity <- function(input_data, model = eADAGEmodel, HW_cutoff = 2.5) {

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
  value_only <- input_data[, -1]

  HWactivity_perGene_pos <- matrix(, nrow = ncol(value_only), ncol = model_size)
  HWactivity_perGene_neg <- matrix(, nrow = ncol(value_only), ncol = model_size)

  for (node in 1:model_size) {

    # positive side
    pos_cutoff <- mean(weight_matrix[, node]) +
                  HW_cutoff * sd(weight_matrix[, node])
    HWG_index <- weight_matrix[, node] >= pos_cutoff
    HWactivity_perGene_pos[, node] <- t(t(weight_matrix[HWG_index, node]) %*%
                                          as.matrix(value_only[HWG_index, ]))
    if (sum(HWG_index) > 0) {
      HWactivity_perGene_pos[, node] <- HWactivity_perGene_pos[, node] /
        sum(HWG_index)
    }

    # negative side
    neg_cutoff <- mean(weight_matrix[, node]) -
                  HW_cutoff * sd(weight_matrix[, node])
    HWG_index <- weight_matrix[, node] <= neg_cutoff
    HWactivity_perGene_neg[, node] <- t(t(weight_matrix[HWG_index, node]) %*%
                                          as.matrix(value_only[HWG_index, ]))
    if (sum(HWG_index) > 0) {
      HWactivity_perGene_neg[, node] <- HWactivity_perGene_neg[, node] /
        sum(HWG_index)
    }
  }

  # combine positive and negative sides
  HWactivity_perGene <- cbind(HWactivity_perGene_pos, HWactivity_perGene_neg)
  rownames(HWactivity_perGene) <- colnames(value_only)
  colnames(HWactivity_perGene) <- c(paste0("Node", seq(1, model_size), "pos"),
                                    paste0("Node", seq(1, model_size), "neg"))

  return(HWactivity_perGene)

}
