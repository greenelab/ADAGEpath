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
#' @export
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


#' Gene set enrichment test
#'
#' Performs a one-side fisher exact test to determine whether two sets of genes
#' have significant overlap.
#'
#' @param set1 character vector storing genes in the first gene set
#' @param set2 character vector storing genes in the second gene set
#' @param set_all character vector storing all possible genes
#' @return pvalue in the fisher exact test
enrich_test <- function(set1, set2, set_all){

  set_overlap <- length(intersect(set1, set2))
  set1_only <- length(set1) - set_overlap
  set2_only <- length(set2) - set_overlap
  others <- length(set_all) - set_overlap - set2_only - set1_only
  contingency_table <- matrix(c(set_overlap, set1_only, set2_only, others),
                              nrow = 2)
  pvalue <- fisher.test(contingency_table, alternative='greater')$p.value

  return(pvalue)
}


#' Signature overlap calculation
#'
#' Tests how significant all possible pairs of signatures in the input list
#' overlap with each other.
#'
#' @param signature_list a list of signatures extracted from an ADAGE model.
#' @return a named list storing adjusted p values for all signature overlap.
#' tests
#' @export
test_signature_overlap <- function(signature_list){

  # generate all possible combinations
  comb_sig <- combn(signature_list, 2)
  comb_name <- combn(names(signature_list), 2)

  # get unique genes in all signatures
  all_HWGs <- unique(unlist(signature_list))

  # test the significance of overlap between every two signature combinations
  # TODO: use parellel version of mapply
  pvalue_list <- mapply(enrich_test, comb_sig[1, ], comb_sig[2, ],
                        MoreArgs = list(set_all = all_HWGs))

  # multiple hypothesis correction
  qvalue_list <- p.adjust(pvalue_list, method = "fdr")
  # name each element with its signature combination
  names(qvalue_list) <- sapply(1: ncol(comb_name),function(x)
    paste(comb_name[, x], collapse = "_"))

  return(qvalue_list)

}


#' Signature overlap matrix construction
#'
#' Converts the list of adjusted p values in the signature overlap tests into
#' a symmetric matrix.
#'
#' @param signature_list a list of signatures extracted from an ADAGE model.
#' @param overlap_qvalues a list of adjusted p values obtained from the function
#' test_signature_overlap()
#' @return a symmetric matrix storing the adjusted p values in the signature
#' overlap tests between every possible signature pair. The diagnol is set to
#' 0 because a signature always has perfect overlap with itself.
#' @export
build_signature_overlap_matrix <- function(signature_list, overlap_qvalues){

  overlap_matrix = matrix(, nrow = length(signature_list),
                          ncol = length(signature_list))
  overlap_matrix[lower.tri(overlap_matrix, diag = FALSE)] <- overlap_qvalues
  overlap_matrix[upper.tri(overlap_matrix)] <-
    t(overlap_matrix)[upper.tri(overlap_matrix)]
  diag(overlap_matrix) <- 0
  rownames(overlap_matrix) <- names(signature_list)
  colnames(overlap_matrix) <- names(signature_list)

  return(overlap_matrix)
}


#' Signature overlap plot
#'
#' Plots the overlap significance between signature pairs using heatmap. Values
#' in the heatmap represent -log10(qvalue). A higher value indicates a more
#' significant overlap. The diagnol is set to the highest value in the selected
#' signature pairs for visualization purpose. The plot will be difficult to
#' see if including more than 50 signatures.
#'
#' @param overlap_matrix a symmetric matrix storing adjusted p values between
#' each signature pair in the overlap test
#' @param signatures a character vector specifying which signatures to include
#' in the plot
#' @export
plot_signature_overlap <- function(overlap_matrix, signatures = NULL){

  if (!is.null(signatures)){
    # make sure all input signatures can be found in the overlap matrix
    if (all(signatures %in% rownames(overlap_matrix))){
      overlap_matrix <- overlap_matrix[signatures, signatures]
    } else {
      stop("Given signatures are not found in the overlap matrix!")
    }

  }

  # set the diagnol to the minimum q value other than 0 in the matrix
  diag(overlap_matrix) <- min(overlap_matrix[overlap_matrix!= 0], na.rm = TRUE)
  # convert the q values into -log10 scale
  overlap_matrix <- -log10(overlap_matrix)

  corrplot::corrplot(overlap_matrix, is.corr = FALSE, method = 'square',
                     order = 'hclust')
}
