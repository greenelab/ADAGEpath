#' Marginal activity calculation
#'
#' Calculates the marginal activities between every two combinations of the input
#' signatures. Marginal activity is defined as the activity of a signature
#' after removing genes that it overlaps with another signature.
#'
#' @param input_data a data.frame with gene IDs in the first column and
#' expression values from the second column.
#' @param selected_signatures a vector storing names of selected signatures
#' @param model  an ADAGE model to extract signatures from
#' @return a data.frame storing marginal activities for the input samples.
#' Its rownames is set as "signature1-signature2",
#' indicating a marginal activity for signature1 after
#' removing the effect of signature2. If rowname is "signature1-signature1",
#' then it's the activity of signature1 itself.
#' @export
calculate_marginal_activity <- function(input_data, selected_signatures,
                                        model) {

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
    # replace NA with 0
    marginal_activities[is.na(marginal_activities)] <- 0

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


#' Marginal significance matrix preparation
#'
#' Converts the limma result on margianl activity into a signatureN * signatureN
#' matrix format.
#'
#' @param marginal_limma_result a data.frame that stores the limma result table
#' returned by the build_limma() function when used on marginal activity.
#' It's rownames is in the format of "signature1-signature2".
#' @return a matrix storing the adjusted p values from
#' the activation test when the effect of the column signature is removed from
#' the row signature. Values in the diagonal of the matrix are the activation
#' significance of signatures themselves.
prepare_marginal_matrix <- function(marginal_limma_result){

  # extract the names of the signature pair
  marginal_limma_result$sig1 <- sapply(rownames(marginal_limma_result),
                                       function(x) unlist(strsplit(x, "-"))[1])
  marginal_limma_result$sig2 <- sapply(rownames(marginal_limma_result),
                                       function(x) unlist(strsplit(x, "-"))[2])

  # long to wide conversion
  marginal_matrix_qvalue <- reshape2::dcast(marginal_limma_result, sig1~sig2,
                                            value.var = "adj.P.Val")
  marginal_matrix_qvalue <- tibble::column_to_rownames(marginal_matrix_qvalue,
                                                       var = "sig1")
  marginal_matrix_qvalue <- as.matrix(marginal_matrix_qvalue)

  return(marginal_matrix_qvalue)
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
#' returned by the build_limma() function when used on marginal activity.
#' It's rownames is in the format of "signature1-signature2".
#' @param signature_order a vector of signature names, the order of signatures
#' in this vector will be used to order signatures in the plot. If NULL,
#' signatures will be ordered alphabatically. (default: NULL)
#' @param sig_cutoff a numeric value used as the significance cutoff.
#' Significance values below the cutoff will be crossed out in the plot.
#' (default: 0.05)
#' @export
plot_marginal_activation <- function(marginal_limma_result,
                                     signature_order = NULL,
                                     sig_cutoff = 0.05){

  marginal_matrix_qvalue <- prepare_marginal_matrix(marginal_limma_result)

  if (!is.null(signature_order)) {
    # reorder matrix with provided signature order
    marginal_matrix_qvalue <- marginal_matrix_qvalue[signature_order,
                                                     signature_order]
  }
  # convert to the negative log10 scale
  marginal_matrix_log10qvalue <- -log10(marginal_matrix_qvalue)

  col <- colorRampPalette(c("white", "yellow", "red"))
  corrplot::corrplot(marginal_matrix_log10qvalue, is.corr = FALSE,
                     col = col(100), p.mat = marginal_matrix_qvalue,
                     sig.level = sig_cutoff, pch.cex = 1)

}


#' Redundant signature removal
#'
#' Removes signatures whose activity changes are no longer significant after
#' removing genes overlapped with another signature.
#'
#' @param marginal_limma_result a data.frame that stores the limma result table
#' returned by the build_limma() function when used on marginal activity.
#' It's rownames is in the format of "signature1-signature2".
#' @param sig_cutoff a numeric value used as the significance cutoff.
#' (default: 0.05)
#' @return a vector storing remaining signatures after removing redundant oens.
#' @export
remove_redundant_signatures <- function(marginal_limma_result,
                                        sig_cutoff = 0.05){

  marginal_matrix_qvalue <- prepare_marginal_matrix(marginal_limma_result)

  binary_matrix <- marginal_matrix_qvalue
  binary_matrix[binary_matrix <= sig_cutoff] <- TRUE
  binary_matrix[binary_matrix != TRUE] <- FALSE

  remove_set <- c()
  save_set <- c()
  # covered signature are signatures that can be covered by other signatures.
  covered_sigs <- which(apply(binary_matrix, 1, all) == FALSE)

  for (i in covered_sigs) {
    # major signatures are signatures that can cover signature i.
    major_sigs <- which(binary_matrix[i, ] == FALSE)

    if (all(major_sigs %in% covered_sigs)) {
      # if all the signatures that can cover signature i are also covered by
      # some other signatures, we loop through them to see how we can
      # remove signature i safely. Signature i can be safely removed only if one
      # of its major signatures is not removed.
      i_pvalue <- marginal_matrix_qvalue[i, i]
      remove_flag <- FALSE
      for (j in major_sigs) {
        j_pvalue <- marginal_matrix_qvalue[j, j]
        if (binary_matrix[j, i]  == FALSE) {
          # signature i and signature j cover each other. In this case, we only
          # remove signature i if it is less significant than signature j. In
          # this way, the more significant signature is kept.
          if(i_pvalue > j_pvalue) {
            remove_set <- c(remove_set, i)
            remove_flag <- TRUE
          }
        } else {
          # signature j covers signature i but is covered by some other
          # signatures. In this case, we remove signature i but save signature j.
          remove_set <- c(remove_set, i)
          remove_flag <- TRUE
          save_set <- c(save_set, j)
        }
        if (remove_flag == TRUE) {
          # stop checking next signature once signature i has been removed
          break
        }
      }

    } else {
      remove_set <- c(remove_set, i)
    }
  }
  # remove signatures in the save set from the remove set
  final_remove_set <- setdiff(remove_set, save_set)
  # get the non-redundant signatures
  non_redundant_sigs <- rownames(marginal_matrix_qvalue)[
    setdiff(1:nrow(marginal_matrix_qvalue), final_remove_set)]

  return(non_redundant_sigs)
}
