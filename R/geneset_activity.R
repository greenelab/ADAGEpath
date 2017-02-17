#' Calculating gene set activity inside signatures
#'
#' Calculates activity of signatures only using genes in the provided gene sets
#'
#' @param signature_geneset_df a data.frame with the first column specifying
#' signature names and second column specifying gene set names.
#' @param gene_set_list a named list storing all possible gene sets with each
#' element being a vector of gene IDs. It is used to extract genes IDs given
#' the name of a gene set.
#' @param input_data a data.frame with gene IDs in the first column and
#' expression values from the second column.
#' @param model the ADAGE model to be used for calculating signature activity.
#' @param HW_cutoff number of standard deviations away from mean in a node's
#' weight distribution to be considered as high-weight (default to 2.5).
#' @return a data.frame with the first column being "signature|geneset" names and
#' the rest columns storing activities for every sample in the input_data.
#' @export
signature_geneset_activity <- function(signature_geneset_df, gene_set_list,
                                       model, input_data, HW_cutoff = 2.5){

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
  express_matrix <- as.matrix(input_data[, -1])

  # calculate activities for each signature-geneset combination
  activities <- sapply(1:nrow(signature_geneset_df), function(x) {
    node <- as.numeric(gsub("\\D", "", signature_geneset_df[x, 1]))
    side <- ifelse(grepl("pos", signature_geneset_df[x, 1]), "pos", "neg")
    gene_set <- gene_set_list[[as.character(signature_geneset_df[x, 2])]]

    one_signature_activity(weight_matrix = weight_matrix,
                           express_matrix = express_matrix,
                           node = node, side = side, gene_set = gene_set,
                           HW_cutoff = HW_cutoff)
  })

  # transpose to have signatures in rows
  activities <- dplyr::as_data_frame(t(activities))

  # omit positive and negative signs of activities
  activities <- abs(activities)

  colnames(activities) <- colnames(input_data)[-1]

  # add the signature name column in the front
  activities <- data.frame(
    signature = paste0(signature_geneset_df[, 1],
                       "|", signature_geneset_df[, 2]),
    activities, stringsAsFactors = FALSE, check.names = FALSE)

  return(activities)

}


#' Combining geneset outputs
#'
#' Combines the signature geneset association output and the geneset limma
#' output.
#'
#' @param signature_geneset_association a data.frame storing signatures and
#' their significantly enriched gene sets, returned by the function
#' annotate_signatures_with_genesets().
#' @param geneset_limma_result a data.frame that stores the result table
#' returned by limma. It includes logFC, adj.P.Val, and other statistics for
#' each feature. The features are in the format "signature|geneset".
#' @return a data.frame that stores both signature geneset association and
#' genesets' limma result.
#' @export
combine_geneset_outputs <- function(signature_geneset_association,
                                    geneset_limma_result){

  # retrieve gene set and signature names
  geneset_limma_result$signature <- sapply(
    rownames(geneset_limma_result), function(x) unlist(strsplit(x, "\\|"))[[1]])
  geneset_limma_result$geneset <- sapply(
    rownames(geneset_limma_result), function(x) unlist(strsplit(x, "\\|"))[[2]])

  # combine limma result with geneset association result
  geneset_limma_result <- geneset_limma_result %>%
    dplyr::select(signature, geneset, logFC, adj.P.Val)
  signature_geneset_association <- signature_geneset_association %>%
    dplyr::select(-pvalue)
  combined_result <- dplyr::left_join(geneset_limma_result,
                                      signature_geneset_association,
                                      by = c("signature", "geneset"))
  combined_result <- combined_result %>% dplyr::rename(
    `activity difference` = logFC, `gene set enrichment (adj.P.val)` = qvalue)
  combined_result <- combined_result[order(combined_result$adj.P.Val), ]

  return(combined_result)
}
