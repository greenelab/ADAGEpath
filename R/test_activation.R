#' Building a linear regression model for differential analysis
#'
#' Uses limma to build a linear model to test the differential expression or
#' differential activation of either gene features or signature features.
#' Only two-group comparison is implemented for now.
#'
#' @param input_data a data.frame that stores either the signature activities or
#' gene expression values. The first column specifies feature names (genes or
#' signatures).
#' @param phenotypes a factor (or a charactor that can be converted into a
#' factor) with two levels that describes the phenotype of each sample.
#' The first factor level is used as the control group and the second level is
#' used as the treatment group in the limma test.
#' @param use.bonferroni a logical value indicating whether to use the more
#' conservative "bonferroni" method in the p value adjustment.
#' This is recommended when there are too many significant features when
#' using the default "BH" method. (default: FALSE)
#' @return a data.frame that stores the result table returned by limma. It
#' includes logFC, adj.P.Val, and other statistics for each feature. LogFC only
#' makes sense when input is gene expression values. If the
#' input is signature activity, then logFC equals to absolute difference.
#' Features with positive logFC have higher values in the second phenotypes
#' level.
#' @seealso \url{https://bioconductor.org/packages/release/bioc/html/limma.html}
#' @export
build_limma <- function(input_data, phenotypes, use.bonferroni = FALSE){

  if (!check_input(input_data)) {
    stop("The input_data should be a data.frame with first column storing
         feature names in character and the rest columns storing
         expression/activity values for each sample in numeric.")
  }

  # rebuild the phenotypes factor if the input phenotypes is a subset of
  # a factor with more than two levels or if it is a character
  phenotypes <- factor(as.character(phenotypes))

  if (nlevels(phenotypes) > 2){
    stop("This function can only deal with two phenotype levels")
  }

  # preperation for limma analysis
  design <- model.matrix(~phenotypes)
  colnames(design) <- levels(phenotypes)

  # convert the first column to rowname as limma use rowname as feature name
  rownames(input_data) <- NULL
  input_data <- tibble::column_to_rownames(input_data,
                                           var = colnames(input_data)[1])

  # use limma to find differentially active features
  fit <- limma::lmFit(input_data, design)
  fit <- limma::eBayes(fit)

  if (use.bonferroni) {
    limma_result <- limma::topTable(fit, coef = 2, number = nrow(input_data),
                                    adjust.method = "bonferroni", sort.by = "none")
  } else {
    limma_result <- limma::topTable(fit, coef = 2, number = nrow(input_data),
                                    adjust.method = "BH", sort.by = "none")
  }

  return(limma_result)
}


#' Most active signature retrieval
#'
#' Returns signatures that are most differentially active in the limma result.
#' We provide three methods to define "most" active: sort by the difference in
#' activities (method = "diff"), filter by a significance level
#' (method = "pvalue"), and consider
#' both difference and significance using pareto fronts (method = "pareto").
#' This function should be used after build_limma() function.
#'
#' @param limma_result a data.frame that stores the limma result table
#' returned by the build_limma() function. It should have rownames that
#' specifies the name of each feature.
#' @param phenotypes a factor (or a charactor that can be converted into a
#' factor) with two levels that describes the phenotype of
#' each sample, should be the same as the phenotypes provided to the
#' build_limma() function.
#' @param method character, can be "diff", "pvalue", or "pareto"
#' (default: "pareto")
#' @param pheno_group character, can be "both", "level1", or "level2".
#' If "both", signatures active in both phenotypes will be merged, sorted, and
#' considered together. If "level1", only signatures active in the first
#' phenotypes level (signatures with negative logFC) will be returned.
#' Similarly, "level2" will return only signatures active in the second
#' phenotypes level (signatures with positive logFC). (default: "both")
#' @param N_signatures int, number of top signatures to return in the "diff"
#' method (default: 10)
#' @param N_fronts int, number of pareto fronts to return in the "pareto" method
#' (default: 5)
#' @param significance_cutoff numeric, significance cutoff
#' to use in the "pvalue" method (default: 0.05).
#' @return a vector storing active signatures.
#' @export
get_active_signatures <- function(limma_result, method = "pareto",
                                  pheno_group = "both",
                                  N_fronts = 5, N_signatures = 10,
                                  significance_cutoff = 0.05) {

  if (!method %in% c("pvalue", "diff", "pareto")){
    stop("Method not recognized! It should be \"pvalue\", \"FC\", or \"pareto\".")
  }

  # build the test_result data.frame with only difference in mean and the adjusted
  # p value after -log10 transform
  test_result <- tibble::data_frame(signature = rownames(limma_result),
                                    diff = limma_result$logFC,
                                    neglog10qvalue = -log10(limma_result$adj.P.Val))

  if (pheno_group == "level1") {
    test_result <- test_result[test_result$diff <= 0, ]
  } else if (pheno_group == "level2") {
    test_result <- test_result[test_result$diff >= 0, ]
  } else if (!pheno_group == "both") {
    stop("pheno_group can only be \"level1\", \"level2\", or \"both\"!")
  }
  test_result$diff <- abs(test_result$diff)

  if (method == "pvalue") {
    active_signatures <- test_result$signature[
      test_result$neglog10qvalue > -log10(significance_cutoff)]
  } else if (method == "diff") {
    active_signatures <- test_result$signature[order(
      test_result$diff, decreasing = TRUE)][1:N_signatures]
  } else if (method == "pareto") {
    # get signatures in the top N layers of pareto fronts of fold change and
    # p values
    active_signatures <- get_paretofront(input_data = test_result,
                                         N_fronts = N_fronts)
  }

  return(active_signatures)
}


#' Pareto fronts calculation
#'
#' Returns solutions in the top N layers of pareto fronts with the goal to
#' maximize two objectives.
#'
#' @param input_data a data.frame with the first column storing a
#' solution's name and second and third columns storing values of its two
#' objectives.
#' @param N_fronts int, number of pareto fronts
#' @return a character storing solutions made into the top N_fronts.
get_paretofront <- function(input_data, N_fronts) {

  # convert to data.frame in case the input_data is a tibble
  input_data <- as.data.frame(input_data)

  # save overlapping solutions into duplicates and then remove duplicates
  duplicates <- input_data[duplicated(input_data[, -1]), ]
  no_dup_data <- input_data[!duplicated(input_data[, -1]), ]

  # order by the second objective and then the first objective
  no_dup_data <- no_dup_data[order(no_dup_data[, 2], no_dup_data[, 3],
                                   decreasing = TRUE), ]

  # get all solutions in the first N pareto fronts
  all_front_sols <- c()
  for (i in 1:N_fronts) {
    # find solutions in the pareto front
    front_sols <- no_dup_data[which(!duplicated(cummax(no_dup_data[, 3]))), 1]
    # remove these solutions from the next loop
    no_dup_data <- dplyr::filter(no_dup_data,
                                 !(no_dup_data[, 1] %in% front_sols))
    all_front_sols <- c(all_front_sols, front_sols)

  }

  # if there are duplicates, check whether dupliates overlap with solutions
  # in pareto fronts
  if (nrow(duplicates) > 0) {

    # combine solutions in pareto fronts and duplicates
    combined_data <- rbind(dplyr::filter(input_data,
                                         input_data[, 1] %in% all_front_sols),
                           duplicates)
    # append solutions in duplicates that overlap with solutions in
    # pareto fronts
    all_front_sols <- c(all_front_sols,
                        combined_data[duplicated(combined_data[, -1]), 1])
  }

  return(all_front_sols)
}


#' Signature volcano plot
#'
#' Makes a volcano plot that shows the activity difference and significance of
#' each signature after testing them with build_limma().
#'
#' @param limma_result a data.frame that stores the limma result table
#' returned by the build_limma() function. It should have rownames that
#' specifies the name of each feature.
#' @param highlight_signatures a character, if provided, signatures in it will
#' be labeled and colored in red in the volcano plot (default: NULL).
#' @param interactive logical, whether the volcano plot should be interactive.
#' If TRUE, the plot is made using plotly. (default: FALSE)
#' @export
plot_volcano <- function(limma_result, highlight_signatures = NULL,
                         interactive = FALSE){

  # build the test_result data.frame with only difference in mean and the adjusted
  # p value after -log10 transform
  test_result <- tibble::data_frame(signature = rownames(limma_result),
                                    diff = limma_result$logFC,
                                    neglog10qvalue = -log10(limma_result$adj.P.Val))

  if (interactive) {
    # significance cutoff, used to plot a dotted line in the volcano plot
    test_result$cutoff <- -log10(0.05)

    if (!is.null(highlight_signatures)) {
      # build a new column to indicate whether a signature belongs to the
      # highlighted signatures. Used to assign colors in the volcano plot.
      test_result$highlight <- ifelse(
        test_result$signature %in% highlight_signatures,"active", "other")
      pal <- c("red", "grey")

      # extract results of the highlighted signatures and only annotate them
      # in the volcano plot
      active_test_result <- test_result[
        test_result$signature %in% highlight_signatures, ]
      highlight_annotation <- list(x = active_test_result$diff,
                                   y = active_test_result$neglog10qvalue,
                                   text = active_test_result$signature,
                                   xref = "x", yref = "y", arrowhead = 0,
                                   ax = 10, ay = -20)

      plotly::plot_ly(data = test_result) %>%
        plotly::add_trace(x = ~diff, y = ~neglog10qvalue,
                          text = ~signature, color = ~highlight, colors = pal,
                          mode = "markers") %>%
        plotly::add_trace(x = ~diff, y = ~cutoff, type = "scatter",
                          mode = "lines", name = "significance cutoff",
                          line = list(color = "red", dash = "dash")) %>%
        plotly::layout(annotations = highlight_annotation,
                       xaxis = list(title = "activity difference"),
                       yaxis = list(title = "significance (-log10qvalue)"))
    } else {

      plotly::plot_ly(data = test_result) %>%
        plotly::add_trace(x = ~diff, y = ~neglog10qvalue, name = "signature",
                          text = ~signature, mode = "markers") %>%
        plotly::add_trace(x = ~diff, y = ~cutoff, type = "scatter",
                          mode = "lines", name = "significance cutoff",
                          line = list(color = "red", dash = "dash")) %>%
        plotly::layout(xaxis = list(title = "activity difference"),
                       yaxis = list(title = "significance (-log10qvalue)"))
    }
  } else {

    plot(test_result$diff, test_result$neglog10qvalue, col = "grey",
         xlab = "activity difference", ylab = "significance (-log10qvalue)",)

    if (!is.null(highlight_signatures)) {

      # only highlight the provided signatures in the plot
      active_test_result <- test_result[
        test_result$signature %in% highlight_signatures, ]
      points(active_test_result$diff, active_test_result$neglog10qvalue,
             pch = 20, col = "red")
      text(active_test_result$diff, active_test_result$neglog10qvalue,
           labels = active_test_result$signature,
           cex = 0.4, pos = 1, offset = 0.3)

    } else {

      # label all signatures
      text(test_result$diff, test_result$neglog10qvalue,
           labels = test_result$signature,
           cex = 0.3, pos = 1, offset = 0.3)
    }

    # add a line to indicate 0.05 significance level
    abline(h = -log10(0.05), col = 'red')
  }
}
