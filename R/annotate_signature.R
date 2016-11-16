#' Gene set fetch
#'
#' Fetches GO or KEGG gene sets from the Tribe webserver
#' (url{tribe.greenelab.com}).
#'
#' @param type character, type of the gene sets, either "GO" or "KEGG".
#' @param max_size int, maximum gene set size to be considered as a meaningful
#' gene set.
#' @param min_size int, minimum gene set size to be considered as a meaningful
#' gene set.
#' @return a named list with each element being a gene set
#' @export
fetch_geneset <- function(type = "KEGG", max_size = 100, min_size = 5){

  # make sure type is one of "GO" or "KEGG"
  if (!type %in% c("GO", "KEGG")) {
    stop("type can only be either GO or KEGG.")
  }

  request_limit <- 2000  # fetch 2000 gene sets at a time
  request_time <- 1
  gene_sets_df_list <- list()
  while (TRUE) {
    request_offset <- (request_time - 1) * request_limit
    tribe_req <- httr::GET("http://tribe.greenelab.com/api/v1/geneset/",
                           query = list(title__startswith = type,
                                        show_tip = "true",
                                        organism = "9",  # P.a. is number 9
                                        xrdb = "Symbol",
                                        limit = request_limit,
                                        offset = request_offset,
                                        format = "json"))

    tribe_content <- suppressWarnings(httr::content(tribe_req, as = "text"))
    gene_sets <- jsonlite::fromJSON(tribe_content)
    gene_sets_df <- dplyr::data_frame(title = gene_sets$objects$title,
                                      count = gene_sets$objects$tip_item_count,
                                      genes = gene_sets$objects$tip$genes)
    gene_sets_df_list[[request_time]] <- gene_sets_df

    # the actual number of returned gene sets
    returned_size <- nrow(gene_sets_df)
    # repeatly request if number of returned gene sets reaches the upper
    # limit of request size, otherwise, break the request loop.
    if (returned_size < request_limit) {
      break
    } else {
      request_time <- request_time + 1
    }
  }
  gene_sets_df <- dplyr::bind_rows(gene_sets_df_list)

  # subset gene sets using max_size and min_size
  gene_sets_df <- gene_sets_df[gene_sets_df$count >= min_size &
                                 gene_sets_df$count <= max_size, ]
  gene_sets_list <- gene_sets_df$genes
  names(gene_sets_list) <- gene_sets_df$title

  return(gene_sets_list)
}


#' Signature annotation
#'
#' Annotates signatures using GO or KEGG terms downloaded from Tribe. Annotation
#' is done through a gene set enrichment test (fisher exact test) with FDR
#' correction.
#'
#' @param signatures a named list of signatures returned by the function
#' extract_signatures().
#' @param annotation_terms a named list of annotation terms returned by the
#' function fetch_geneset().
#' @param all_genes a vector storing all possible background genes in the
#' enrichment test.
#' @param significance_cutoff numeric, FDR significance cutoff used to filter
#' the test result (default to 0.05).
#' @return a data.frame storing significantly enriched annotation terms for all
#' signatures.
#' @export
annotate_signatures <- function(signatures, annotation_terms, all_genes,
                                significance_cutoff = 0.05){

  # enrichment test between every signature - term pair
  pvalue_matrix <- sapply(signatures, function(y)
    sapply(annotation_terms, function(x) enrich_test(x, y, all_genes)))
  colnames(pvalue_matrix) <- names(signatures)
  rownames(pvalue_matrix) <- names(annotation_terms)

  pvalue_melted <- reshape2::melt(t(pvalue_matrix))
  colnames(pvalue_melted) <- c("signature", "term", "pvalue")

  pvalue_melted$qvalue <- p.adjust(pvalue_melted$pvalue, method = "fdr")
  pvalue_melted_sig <- pvalue_melted[pvalue_melted$qvalue <=
                                       significance_cutoff, ]

  pvalue_melted_sig <- pvalue_melted_sig[order(pvalue_melted_sig$signature), ]

  return(pvalue_melted_sig)
}
