#' Signature extraction
#'
#' Extracts all gene signatures from an ADAGE model.
#'
#' @param model the ADAGE model to be used for extracting gene signatures
#' (default: the 300-node eADAGE model preloaded in the package).
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
#' @return value returned by fisher.test() function, which is a list with
#' class "htest" containing p.value, conf.int, estimate, and so on.
#' @seealso \code{\link[stats]{fisher.test}}
#' @export
enrich_test <- function(set1, set2, set_all){

  set_overlap <- length(intersect(set1, set2))
  set1_only <- length(set1) - set_overlap
  set2_only <- length(set2) - set_overlap
  others <- length(set_all) - set_overlap - set2_only - set1_only
  contingency_table <- matrix(c(set_overlap, set1_only, set2_only, others),
                              nrow = 2)
  fisher_result <- fisher.test(contingency_table, alternative='greater')

  return(fisher_result)
}


#' Signature overlap calculation
#'
#' Tests how significant all possible pairs of input signatures
#' overlap with each other in term of their gene compositions.
#'
#' @param selected_signatures a vector storing names of signatures selected
#' to be tested.
#' @param model an ADAGE model to extract signatures from
#' (default: the 300-node eADAGE model preloaded in the package).
#' @return a named list storing odds ratios for all pairs of signature overlap
#' tests.
#' @export
test_signature_overlap <- function(selected_signatures, model = eADAGEmodel){

  # extract all signatures' gene lists
  signatures_genes <- extract_signatures(model)
  # get unique genes in all signatures
  all_genes <- unique(unlist(signatures_genes))
  # make sure the input selected_signatures can be found in the model
  if (!all(selected_signatures %in% names(signatures_genes))){
    stop("Selected signatures can not be found in the specified model.")
  }
  # get the gene lists of only the selected signatures
  selected_signatures_genes <- signatures_genes[selected_signatures]

  # generate all possible combinations of both signatures' gene lists and names
  comb_sig <- combn(selected_signatures_genes, 2)
  comb_name <- combn(names(selected_signatures_genes), 2)

  # test the significance of overlap between every two signature combinations
  result_table <- mapply(enrich_test, comb_sig[1, ], comb_sig[2, ],
                        MoreArgs = list(set_all = all_genes))

  # name each element with its signature combination
  colnames(result_table) <- sapply(1: ncol(comb_name),function(x)
    paste(comb_name[, x], collapse = "_"))

  # extract odds ratio from the enrichment results
  odds_ratios <- result_table["estimate", ]

  return(odds_ratios)

}


#' Signature overlap matrix construction
#'
#' Converts the list of odds ratios in the signature overlap tests into
#' a symmetric matrix.
#'
#' @param selected_signatures a vector storing names of selected signatures.
#' @param odds_ratios a list of odds ratios obtained from the function
#' test_signature_overlap()
#' @return a symmetric matrix storing the odds ratios in the signature
#' overlap tests between every possible signature pair. The diagonal is set to
#' Inf because a signature always has perfect overlap with itself.
build_signature_overlap_matrix <- function(selected_signatures, odds_ratios){

  # initialize the matrix
  overlap_matrix = matrix(, nrow = length(selected_signatures),
                          ncol = length(selected_signatures))

  # assign odds ratios to the lower triangle
  overlap_matrix[lower.tri(overlap_matrix, diag = FALSE)] <- unlist(odds_ratios)

  # copy the lower triangle to the upper triangle
  overlap_matrix[upper.tri(overlap_matrix)] <-
    t(overlap_matrix)[upper.tri(overlap_matrix)]

  # add row and column names
  rownames(overlap_matrix) <- selected_signatures
  colnames(overlap_matrix) <- selected_signatures

  # set diagnal to infinity
  diag(overlap_matrix) <- Inf

  return(overlap_matrix)
}


#' Signature overlap plot
#'
#' Plots the overlap significance between signature pairs using a heatmap.
#' Values in the heatmap represent the odds ratio in fisher exact test.
#' A higher odds ratio indicates a more significant overlap between genes
#' in two signatures. The diagonal is set to the highest odds ratio in the
#' selected signature pairs just for visualization purpose. The plot will be
#' difficult to see if including more than 50 signatures.
#'
#' @param selected_signatures a vector storing names of signatures to include
#' in the plot
#' @param model an ADAGE model to extract signatures from
#' (default: the 300-node eADAGE model preloaded in the package).
#' @export
plot_signature_overlap <- function(selected_signatures, model = eADAGEmodel){

  # test how sigfinicant each signature pair overlaps
  odds_ratios <- test_signature_overlap(selected_signatures, model)

  # convert the resulting odds ratio list into a symmetric matrix
  overlap_matrix <- build_signature_overlap_matrix(selected_signatures,
                                                   odds_ratios)

  # for visualization purpose, set the Inf to the maximum odds
  # ratio in the matrix
  overlap_matrix[overlap_matrix == Inf] <-
    max(overlap_matrix[overlap_matrix != Inf], na.rm = TRUE)

  # plot the odds ratio heatmap using corrplot
  corrplot::corrplot(overlap_matrix, is.corr = FALSE, method = 'square',
                     order = 'hclust')
}
