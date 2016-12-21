#' Marginal activity calculation
#'
#' Calculates the marginal activities between all possible pairs of the input
#' signatures. Marginal activity is defined as the activity of a signature
#' after removing genes that overlap with another signature.
#'
#' @param input_data a data.frame with gene IDs in the first column and
#' expression values from the second column.
#' @param selected_signatures a vector storing names of selected signatures
#' @param model  an ADAGE model to extract signatures from
#' (default: the 300-node eADAGE model preloaded in the package).
#' @return a data.frame storing marginal activities for the input samples.
#' It's rownames is set as "signature1-signature2",
#' indicating a marginal activity for signature1 after
#' removing the effect of signature2. If rowname is "signature1-signature1",
#' then it's the activity of signature1.
#' @export
calculate_marginal_activity <- function(input_data, selected_signatures,
                                        model = eADAGEmodel) {

  if (!check_input(input_data)){
    stop("The input data should be a data.frame with first column storing
         geneIDs in character and the rest columns storing expression values
         for each sample in numeric.")
  }

  if (!check_input(model)) {
    stop("The model should be a data.frame with first column being gene IDs
         in character and the rest columns storing numeric weight values for
         each node per column.")
  }

  if(!all(input_data[, 1] == model[, 1])){
    stop("The gene IDs stored in the first column of the input data and the
         first column of the model should be the same.")
  }

  # prepare model's weight matrix and the input data
  weight_matrix <- as.matrix(model[, -1])
  rownames(weight_matrix) <- model[[1]]
  value_only <- as.matrix(input_data[, -1])
  rownames(value_only) <- input_data[[1]]

  # extract the selected signatures from the model
  signatures_genes <- extract_signatures(model)
  selected_signatures_genes <- signatures_genes[selected_signatures]

  # create an index to every two combinations of the signatures
  comb_index <- combn(seq_along(selected_signatures_genes), 2)

  # calculate marginal activities
  marginal_activity_list <- mapply(function(x, y){

    # get genes only in one signature but not the other
    sig1_only <- setdiff(selected_signatures_genes[[x]],
                         selected_signatures_genes[[y]])
    sig2_only <- setdiff(selected_signatures_genes[[y]],
                         selected_signatures_genes[[x]])

    # extract signature index from its name, e.x. 53 from "Node53pos"
    sig1_index <- as.numeric(unique(unlist(regmatches(
      names(selected_signatures_genes)[x],
      gregexpr("[0-9]+",names(selected_signatures_genes)[x])))))
    sig2_index <- as.numeric(unique(unlist(regmatches(
      names(selected_signatures_genes)[y],
      gregexpr("[0-9]+", names(selected_signatures_genes)[y])))))

    # calculate marginal activity using genes in signature only
    sig1_activities <- t(weight_matrix[sig1_only, sig1_index]) %*%
      value_only[sig1_only, ] / length(sig1_only)
    sig2_activities <- t(weight_matrix[sig2_only, sig2_index]) %*%
      value_only[sig2_only, ] / length(sig2_only)

    # name the activity as "sig1-sig2", indicating the activity is calculated
    # using genes exclusively in sig1 but not sig2
    rownames(sig1_activities) <- paste(names(selected_signatures_genes)[c(x,y)],
                                       collapse = "-")
    rownames(sig2_activities) <- paste(names(selected_signatures_genes)[c(y,x)],
                                       collapse = "-")

    marginal_activities <- data.frame(rbind(sig1_activities, sig2_activities),
                                      check.names = FALSE)
    return(marginal_activities)
  }, comb_index[1, ], comb_index[2, ], SIMPLIFY = FALSE)

  marginal_activities <- do.call(rbind, marginal_activity_list)
  marginal_activities <- tibble::rownames_to_column(marginal_activities,
                                                    var = "signature")

  # get the activity of the selected signatures
  full_activities <- calculate_activity(input_data, model)
  selected_activities <- full_activities[match(selected_signatures,
                                               full_activities$signature),]
  # rename signature as "signature-signature" for consistence
  selected_activities$signature <- sapply(selected_activities$signature,
                                          function(x) paste(x, x, sep = "-"))
  # combine with the marginal activities
  combined_activities <- rbind(selected_activities, marginal_activities)

  return(combined_activities)
}


#' Marginal activation plot
#'
#' Plots the activation significance of the marginal effects of signatures.
#' The value of the heatmap represents -log10 transformed adjusted p value from
#' the activation test when the effect of the column signature is removed from
#' the row signature. Values in the diagonal of the heatmap are the activation
#' significance of signatures themselves.
#'
#' @param marginal_limma_result a data.frame that stores the limma result table
#' returned by the build_limma() function. It's rownames are in the format of
#' "signature1-signature2".
#' @export
plot_marginal_activation <- function(marginal_limma_result){

  # extract the names of the signature pair
  marginal_limma_result$sig1 <- sapply(rownames(marginal_limma_result),
                                       function(x) unlist(strsplit(x, "-"))[1])
  marginal_limma_result$sig2 <- sapply(rownames(marginal_limma_result),
                                       function(x) unlist(strsplit(x, "-"))[2])
  marginal_limma_result$log10qvalue <- -log10(marginal_limma_result$adj.P.Val)

  # long to wide conversion
  marginal_matrix <- reshape2::dcast(marginal_limma_result, sig1~sig2,
                                     value.var = "log10qvalue")
  marginal_matrix <- tibble::column_to_rownames(marginal_matrix, var = "sig1")

  col <- colorRampPalette(c("white", "yellow","red"))
  corrplot::corrplot(as.matrix(marginal_matrix), is.corr = FALSE, order = "FPC",
                     col = col(100))

}

