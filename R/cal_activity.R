#' Signature activity calculation
#'
#' Calculates activities for each signature in an ADAGE model
#' for each sample in the input_data
#'
#' @param input_data a data.frame with gene IDs in the first column and
#' expression values from the second column.
#' @param model the ADAGE model to be used for calculating signature activity
#' (default: the 300-node eADAGE model preloaded in the package).
#' @param HW_cutoff number of standard deviations away from mean in a node's
#' weight distribution to be considered as high-weight (default to 2.5).
#' Signature activities are calculated only using HW genes.
#' @return a data.frame with the first column being signature names and
#' the rest columns storing signature activities for every sample in the
#' input_data.
#' @export
calculate_activity <- function(input_data, model = eADAGEmodel,
                               HW_cutoff = 2.5){

  if (!check_input(input_data)){
    stop("The input data should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric
         expression values for each sample.")
  }

  if (!check_input(model)) {
    stop("The model should be a data.frame with the first column as a character
         of gene IDs and the rest of the columns storing numeric weight values
         for each node.")
  }

  if(!all(input_data[, 1] == model[, 1])){
    stop("The gene IDs stored in the first column of the input data and the
         first column of the model should be the same.")
  }

  weight_matrix <- as.matrix(model[, -1])
  rownames(weight_matrix) <- model[[1]]
  model_size <- ncol(weight_matrix)
  express_matrix <- as.matrix(input_data[, -1])

  HWactivity_perGene_pos <- sapply(1:model_size, function(x)
    one_signature_activity(weight_matrix = weight_matrix,
                           express_matrix = express_matrix,
                           node = x, side = "pos", HW_cutoff = HW_cutoff))
  HWactivity_perGene_neg <- sapply(1:model_size, function(x)
    one_signature_activity(weight_matrix = weight_matrix,
                           express_matrix = express_matrix,
                           node = x, side = "neg", HW_cutoff = HW_cutoff))

  # combine positive and negative sides
  HWactivity_perGene <- dplyr::bind_cols(data.frame(HWactivity_perGene_pos),
                                         data.frame(HWactivity_perGene_neg))

  # transpose to have signatures in rows
  HWactivity_perGene <- dplyr::as_data_frame(t(HWactivity_perGene))

  # omit positive and negative signs of activities
  HWactivity_perGene <- abs(HWactivity_perGene)

  colnames(HWactivity_perGene) <- colnames(express_matrix)

  # add the signature name column in the front
  HWactivity_perGene <- data.frame(
    signature = c(paste0("Node", seq(1, model_size), "pos"),
                  paste0("Node", seq(1, model_size), "neg")),
    HWactivity_perGene, stringsAsFactors = FALSE, check.names = FALSE)

  return(HWactivity_perGene)
}


#' One signature activity
#'
#' Calculates activities for a specific signature in an ADAGE model. If gene_set
#' is provided, it will only calculate the acitivity of a signature using
#' genes in the gene_set.
#'
#' @param weight_matrix a data.matrix storing the weight matrix in an ADAGE model.
#' @param express_matrix a data.matrix storing gene expression values in a dataset
#' @param node an int ranging from 1 to number of columns in weight_matrix
#' @param side character, "pos" or "neg"
#' @param gene_set a character vector storing gene IDs, must match
#' the gene IDs used in the model. (default: NULL)
#' @param HW_cutoff number of standard deviations from mean in a node's weight
#' distribution to be considered as high-weight (default to 2.5).
#' @return a vector storing activities of one signature across samples
#' @export
one_signature_activity <- function(weight_matrix, express_matrix, node, side,
                                   gene_set = NULL, HW_cutoff = 2.5){

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

  if (is.null(gene_set)) {
    GS_index <- rep(TRUE, nrow(weight_matrix))
  } else {
    # only include genes in the gene set
    GS_index <- rep(FALSE, nrow(weight_matrix))
    GS_index[match(gene_set, rownames(weight_matrix))] <- TRUE
  }

  combined_index <- HWG_index & GS_index

  node_activity <- t(weight_matrix[combined_index, node]) %*%
    express_matrix[combined_index, ]

  # only divide by number of genes when it is higher than zero
  if (any(combined_index) > 0) {
    node_activity <- node_activity / sum(combined_index)
  }

  return(node_activity)
}


#' Signature activity heatmap
#'
#' Plots a heatmap showing signature activities in a dataset.
#'
#' @param activity a data.frame that stores the signature activities for each
#' sample in a dataset. The first column is signature name and activity values
#' start from the second column.
#' @param signatures a character vector specifying which signatures to include
#' in the heatmap (default: NULL, all signatures will be included).
#' @param fix_color_range logical. If TRUE, fix the heatmap color
#' range to the maximum activity range of the input activity. If FALSE, heatmap
#' color range is determined by the activity ranges of the selected signatures.
#' (default: TRUE)
#' @export
plot_activity_heatmap <- function(activity, signatures = NULL,
                                  fix_color_range = TRUE){

  if (!check_input(activity)){
    stop("The input activity should be a data.frame with the first column as
         a character of signature names and the rest of the columns storing
         numeric activity values for each sample.")
  }

  # set the activity rownames to NULL in case it has rownames
  rownames(activity) <- NULL
  # convert the signature name column to rowname to make it easy to plot heatmap
  activity <- tibble::column_to_rownames(activity, var = "signature")

  # keep a record of maximum activity range of all signatures
  max.activity.range <- max(apply(activity, 1, function(x) diff(range(x))))

  # subset the activity data.frame to only contain the input signatures
  if (!is.null(signatures)) {
    # make sure all input signatures can be found in the overlap matrix
    if (all(signatures %in% rownames(activity))) {
      activity <- activity[signatures, ]
    } else {
      stop("Given signatures are not found in the rownames of the
           activity data.frame!")
    }
  }

  # get the minimum value of each signature activity
  min_val <- apply(activity, 1, min)

  # transform each signature's activity values by subtracting its minimum
  activity_scaled <- t(scale(t(activity), center = min_val, scale = FALSE))

  # activity heatmap color panel
  activity.color <- gplots::colorpanel(255, rgb(0/255, 176/255, 240/255),
                                       rgb(230/255, 230/255, 230/255),
                                       rgb(255/255, 255/255, 0/255))
  if (fix_color_range) {
    # use the maximum activity range of all signatures
    color.range <- seq(0, max.activity.range, length = 256)
  } else {
    # use the maximum activity range of the selected signatures
    color.range <- seq(0, max(activity_scaled), length = 256)
  }

  # plot activity heatmap
  suppressWarnings(gplots::heatmap.2(activity_scaled, Rowv = TRUE, Colv = FALSE,
                    trace = "none", margins = c(10, 20), cexRow = 1, cexCol = 1,
                    col = activity.color, breaks = color.range))
}
