#' Defining sample phenotypes
#'
#' Defines the pehnotype of each sample in a dataset.
#'
#' @param phenotypes a character with phenotypes separated by comma,
#' e.g. "wt,wt,wt,mt,mt,mt"
#' @return a factor storing the input phenotypes
#' @export
set_phenotype <- function(phenotypes){

  phenotypes <- factor(unlist(strsplit(phenotypes, ",")))
  return(phenotypes)
}

#' Linear model for activation detection
#'
#' Uses limma to build a linear model of signature activities and test their
#' activations. Only two-group comparison is implemented for now.
#'
#' @param activity a data.frame that stores the signature activities for each
#' sample in a dataset.
#' @param phenotypes a factor obtained from the set_phenotype() function.
#' @param sample_indices a int vector to specify which samples to include
#' in the linear model.
#' @return a list including test_result and phenotypes: test_result is a
#' data.frame that stores the activation differences and significance;
#' phenotypes is a factor after subseting the input phenotypes with
#' the sample_indices.
#' @export
build_limma <- function(activity, phenotypes, sample_indices = NULL){

  # subset the input activity with specified sample indices
  if (!is.null(sample_indices)) {
    activity <- activity[, sample_indices]
    phenotypes <- phenotypes[sample_indices]
  }

  if (nlevels(phenotypes) > 2){
    stop("This function can only deal with two phenotype groups at this moment.")
  }

  # preperation for limma analysis
  design <- model.matrix(~phenotypes)
  colnames(design) <- levels(phenotypes)

  # use limma to find differentially active signatures
  fit <- limma::lmFit(activity, design)
  fit <- limma::eBayes(fit)
  result <- limma::topTable(fit, coef = 2, number = nrow(activity),
                            adjust.method = "BH", sort.by = "none")

  # use bonferroni method to do the multiple hypothesis correction
  result$Bon.p <- p.adjust(result$P.Value, method = "bonferroni")

  # build the result data.frame with only difference in mean and the adjusted
  # p value after -log10 transform
  sig_FC_pvalue <- tibble::data_frame(signature = rownames(result),
                                      diff = result$logFC,
                                      neglog10qvalue = -log10(result$Bon.p))

  return(list(test_result = sig_FC_pvalue, phenotypes = phenotypes))
}


#' Most active signature retrieval
#'
#' Returns the signatures that are most active in the limma test result. We
#' provide three methods to define "most" active: sort by the activity
#' difference (diff), significance (pvalue), and consider both difference and
#' significance using pareto fronts (pareto).
#'
#' @param limma_result list, returned by the build_limma() function.
#' @param method character, can be "diff", "pvalue", or "pareto"(default)
#' @param pheno_group character, can be "both" (default) or one of the
#' phenotype level used during build_limma(). If "both", signatures active
#' in both phenotypes will be merged, sorted, and returned together. Otherwise,
#' only signatures active in the input phenotype will be returned.
#' @param N_signatures int, number of top signatures to return in the "diff"
#' method (default to 10)
#' @param N_fronts int, number of pareto fronts to return in the "pareto" method
#' (default to 5)
#' @param significance_cutoff numeric, significance cutoff (in the -log10 scale)
#' to use in the "pvalue" method
#' @return a vector storing active signatures
#' @export
get_active_signatures <- function(limma_result, method = "pareto",
                                  pheno_group = "both", N_signatures = 10,
                                  N_fronts = 5,
                                  significance_cutoff = -log10(0.05)) {

  if (!method %in% c("pvalue", "diff", "pareto")){
    stop("Method not recognized! It should be \"pvalue\", \"FC\", or \"pareto\".")
  }

  test_result <- limma_result$test_result
  phenotypes <- limma_result$phenotypes

  if (pheno_group == levels(phenotypes)[1]) {
    test_result <- test_result[test_result$diff <= 0, ]
  } else if (pheno_group == levels(phenotypes)[2]) {
    test_result <- test_result[test_result$diff >= 0, ]
  } else if (!pheno_group == "both") {
    stop("pheno_group can only be \"both\" or one of the two phenotypes used
         in the run_limma function!")
  }
  test_result$diff <- abs(test_result$diff)

  if (method == "pvalue") {
    active_signatures <- test_result$signature[
      test_result$neglog10qvalue > significance_cutoff]
  } else if (method == "diff") {
    active_signatures <- test_result$signature[order(
      test_result$diff, decreasing = TRUE)][1:N_signatures]
  } else if (method == "pareto") {
    # get signatures in the top N layers of pareto fronts of fold change and
    # p values
    active_signatures <- get_paretofront(test_result = test_result,
                                         N_fronts = N_fronts)
  }

  return(active_signatures)

}


#' Pareto fronts calculation
#'
#' @param test_result a data.frame included in the list returned from
#' the build_limma() function.
#' @param N_fronts int, number of pareto fronts
#' @return a charactor storing all signatures made into the top N_fronts.
get_paretofront <- function(test_result, N_fronts) {

  # convert to data.frame in case the test_result is a tibble
  test_result <- as.data.frame(test_result)

  # save overlapping signatures into duplicates and then remove duplicates
  duplicates <- test_result[duplicated(test_result[, -1]), ]
  no_dup_data <- test_result[!duplicated(test_result[, -1]), ]

  # order by the second objective and then the first objective
  no_dup_data <- no_dup_data[order(no_dup_data[, 2], no_dup_data[, 3],
                                   decreasing = TRUE), ]

  # get all signatures in the first N pareto fronts
  all_front_sigs <- c()
  for (i in 1:N_fronts) {
    # find signatures in the pareto front
    front_sigs <- no_dup_data[which(!duplicated(cummax(no_dup_data[, 3]))), 1]
    # remove these signatures from the next loop
    no_dup_data <- dplyr::filter(no_dup_data,
                                 !(no_dup_data[, 1] %in% front_sigs))
    all_front_sigs <- c(all_front_sigs, front_sigs)

  }

  # if there are duplicates, check whether dupliates overlap with signatures
  # in pareto fronts
  if (nrow(duplicates) > 0) {

    # combine signatures in pareto fronts and duplicates
    combined_data <- rbind(dplyr::filter(test_result,
                                         test_result[, 1] %in% all_front_sigs),
                           duplicates)
    # append signatures in duplicates that overlap with signatures in
    # pareto fronts
    all_front_sigs <- c(all_front_sigs,
                        combined_data[duplicated(combined_data[, -1]), 1])
  }

  return(all_front_sigs)
}


#' Signature volcano plot
#'
#' Makes a volcano plot showing the activity difference and significance of
#' each signature after testing them with limma.
#'
#' @param limma_result list, returned by the build_limma() function.
#' @param active_signatures a character, if provided, signatures in it will
#' be labeled and colored in red in the plot (default to NULL).
#' @export
plot_volcano <- function(limma_result, active_signatures = NULL){

  test_result <- limma_result$test_result
  plot(test_result$diff, test_result$neglog10qvalue, col = "grey",
       xlab = "activity difference", ylab = "significance (-log10qvalue)")

  if (!is.null(active_signatures)) {
    active_test_result <- test_result[
      test_result$signature %in% active_signatures, ]
    points(active_test_result$diff, active_test_result$neglog10qvalue,
           pch = 20, col = "red")
    text(active_test_result$diff, active_test_result$neglog10qvalue,
         labels = active_test_result$signature,
         cex = 0.4, pos = 1, offset = 0.3)
  } else {
    text(test_result$diff, test_result$neglog10qvalue,
         labels = test_result$signature,
         cex = 0.3, pos = 1, offset = 0.3)
  }
  abline(h = -log10(0.05), col = 'red')
}
