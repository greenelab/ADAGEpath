
run_limma <- function(activity, sample_indices = NULL, pheno_groups){
  # pre-process the file to only contain 2 experimental conditions
  if (!is.null(sample_indices)) {
    activity <- activity[, sample_indices]
  }

  # preperation for limma anaysis
  pheno_groups <- factor(unlist(strsplit(pheno_groups, ",")))
  design <- model.matrix(~pheno_groups)
  colnames(design) <- levels(pheno_groups)

  # use limma to find differentially active signatures between two conditions
  fit <- limma::lmFit(activity, design)
  fit <- limma::eBayes(fit)
  # TODO: other function than topTable
  result <- limma::topTable(fit, coef = 2, number = nrow(activity),
                            adjust.method = "BH", sort.by = "none")
  result$Bon.p <- p.adjust(result$P.Value, method = "bonferroni")
  # build a data frame for pareto front analysis on fold change and p value
  sig_FC_pvalue <- tibble::data_frame(signature = rownames(result),
                                      logFC = result$logFC,
                                      neglog10qvalue = -log10(result$Bon.p))
  return(sig_FC_pvalue)
}


get_active_signatures <- function(limma_result, method = "pareto",
                                  two_side = TRUE, top_signatures = 10,
                                  pareto_layers = 5,
                                  neglog10qvalue_cutoff = -log10(0.05)) {

  if (!method %in% c("pvalue", "FC", "pareto")){
    stop("Method not recognized! It should be \"pvalue\", \"FC\", or \"pareto\".")
  }

  if (two_side) {
    limma_result$neglog10qvalue <- abs(limma_result$neglog10qvalue)
  }

  if (method == "pvalue") {
    active_signatures <- limma_result$signature[
      limma_result$neglog10qvalue > neglog10qvalue_cutoff]
  } else if (method == "FC") {
    active_signatures <- limma_result$signature[order(
      limma_result$logFC, decreasing = TRUE)][1:top_signatures]
  } else if (method == "pareto") {
    # get signatures in the top N layers of pareto fronts of fold change and
    # p values
    active_signatures <- get_paretofront(input_data = limma_result,
                                         N_layers = pareto_layers)
  }

  return(active_signatures)

}


get_paretofront <- function(input_data, N_layers) {

  # convert to data.frame in case the input_data is a tibble
  input_data <- as.data.frame(input_data)

  # save overlapping signatures into duplicates and then remove duplicates
  duplicates <- input_data[duplicated(input_data[, -1]), ]
  no_dup_data <- input_data[!duplicated(input_data[, -1]), ]

  # order by the second objective and then the first objective
  no_dup_data <- no_dup_data[order(no_dup_data[, 2], no_dup_data[, 3],
                                   decreasing = TRUE), ]

  # get all signatures in the first N pareto fronts
  all_front_sigs <- c()
  for (i in 1:N_layers) {
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
    combined_data <- rbind(dplyr::filter(input_data,
                                         input_data[, 1] %in% all_front_sigs),
                           duplicates)
    # append signatures in duplicates that overlap with signatures in
    # pareto fronts
    all_front_sigs <- c(all_front_sigs,
                        combined_data[duplicated(combined_data[, -1]), 1])
  }

  return(all_front_sigs)
}


plot_volcano <- function(limma_result, active_signatures = NULL){
  plot(limma_result$logFC, limma_result$neglog10qvalue, type = "n")
  text(limma_result$logFC, limma_result$neglog10qvalue,
       labels = limma_result$signature, cex = 0.4, pos = 1, offset = 0.3)
  if (!is.null(active_signatures)) {
    points(limma_result$logFC[limma_result$signature %in% active_signatures],
           limma_result$neglog10qvalue[limma_result$signature %in% active_signatures],
           pch = 20, col = "red")
  }
  abline(h = -log10(0.05), col = 'red')
}
