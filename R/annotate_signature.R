fetch_geneset <- function(type = "KEGG", max_size = 100, min_size = 5){

  # make sure type is one of "GO" or "KEGG"
  if (!type %in% c("GO", "KEGG")) {
    stop("type can only be either GO or KEGG.")
  }

  request_limit <- 2000
  request_time <- 1
  gene_sets_df_list <- list()
  while (TRUE) {
    request_offset <- (request_time - 1) * request_limit
    tribe_req <- httr::GET("http://tribe.greenelab.com/api/v1/geneset/",
                           query = list(title__startswith = type,
                                        show_tip = "true",
                                        organism = "9",
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
    returned_size <- nrow(gene_sets_df)
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

  return(gene_sets_df)
}


annotate_signatures <- function(signatures, annotation_terms, all_genes,
                                significance_cutoff = 0.05){

  pvalue_matrix <- sapply(signatures, function(y)
    sapply(annotation_terms$genes, function(x) enrich_test(x, y, all_genes)))
  colnames(pvalue_matrix) <- names(signatures)
  rownames(pvalue_matrix) <- annotation_terms$title
  pvalue_melted <- reshape2::melt(t(pvalue_matrix))
  colnames(pvalue_melted) <- c("signature", "term", "pvalue")

  pvalue_melted$qvalue <- p.adjust(pvalue_melted$pvalue, method = "fdr")
  pvalue_melted_sig <- pvalue_melted[pvalue_melted$qvalue <=
                                       significance_cutoff, ]
  pvalue_melted_sig <- pvalue_melted_sig[order(pvalue_melted_sig$signature), ]

  return(pvalue_melted_sig)
}
