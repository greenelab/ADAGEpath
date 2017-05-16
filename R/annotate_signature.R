#' Gene set fetch
#'
#' Fetches gene sets from the Tribe webserver
#' (url{https://tribe.greenelab.com}). It can retrieve GO, KEGG, or user-created
#' public gene sets for Pseudomonas aeruginosa.
#'
#' @param type character, type of the gene sets, can be "GO", "KEGG", or "User"
#' (default: "KEGG").
#' @param username character, creator's TRIBE username (default: NULL)
#' @param access_date character, in the format of month-day-year (12-31-16).
#' When specified, only the version of a term edited right before the access
#' date is retrieved. In default, a term's most up-to-date version is retrieved.
#' (default: NULL)
#' @param max_size int, maximum gene set size to be considered as a meaningful
#' gene set (default: Inf).
#' @param min_size int, minimum gene set size to be considered as a meaningful
#' gene set (default: 0).
#' @return a named list with each element being a gene set
#' @export
fetch_geneset <- function(type = "KEGG", username = NULL, access_date = NULL,
                          max_size = Inf, min_size = 0){

  # make sure type is one of "GO", "KEGG", or "User"
  if (!type %in% c("GO", "KEGG", "User")) {
    stop("type can only be GO or KEGG or User.")
  }

  # username cannot be NULL if type is User
  if (type == "User" & is.null(username)) {
    stop("Please provide the TRIBE username.")
  }

  # to prevent overloading the TRIBE webserver, we limit the download size
  # to 200 gene sets at one time and make multiple downloads through changing
  # the offset until all gene sets have been downloaded.
  request_limit <- 200
  request_time <- 1
  gene_sets_df_list <- list()

  while (TRUE) {
    # set request offset to the number of request limits that have already
    # been reached
    request_offset <- (request_time - 1) * request_limit

    if (is.null(access_date)) {
      # if access_date is NULL, use the tip version

      if (type == "User") {
        # filter gene set by creator's username
        tribe_req <- httr::GET("https://tribe.greenelab.com/api/v1/geneset/",
                               query = list(creator__username = username,
                                            show_tip = "true",
                                            organism = "9",  # P.a. is number 9
                                            xrdb = "Symbol",
                                            limit = request_limit,
                                            offset = request_offset,
                                            format = "json"))
      } else {
        # filter gene set by title
        tribe_req <- httr::GET("https://tribe.greenelab.com/api/v1/geneset/",
                               query = list(title__startswith = type,
                                            show_tip = "true",
                                            organism = "9",  # P.a. is number 9
                                            xrdb = "Symbol",
                                            limit = request_limit,
                                            offset = request_offset,
                                            format = "json"))
      }

      tribe_content <- suppressMessages(httr::content(tribe_req, as = "text"))
      gene_sets <- jsonlite::fromJSON(tribe_content)

      if (gene_sets$meta$total_count > 0) {
        gene_sets_df <- dplyr::data_frame(title = gene_sets$objects$title,
                                          genes = gene_sets$objects$tip$genes,
                                          count = gene_sets$objects$tip_item_count)
        gene_sets_df_list[[request_time]] <- gene_sets_df

        # the actual number of returned gene sets
        returned_size <- nrow(gene_sets_df)
        # If the return_size reaches the request limit, it indicates that we haven't
        # downloaded all gene sets yet. We will make a new request repeatly until
        # the number of returned gene sets is smaller than the request limit.
        if (returned_size < request_limit) {
          break
        } else {
          request_time <- request_time + 1
        }
      } else {
        break
      }
    } else {
      # if access_date is provided, retrieve the gene set version
      # created right before the access_date

      if (type == "User") {
        # filter gene set by creator's username
        tribe_req <- httr::GET("https://tribe.greenelab.com/api/v1/geneset/",
                               query = list(creator__username = username,
                                            show_versions = "true",
                                            organism = "9",  # P.a. is number 9
                                            xrdb = "Symbol",
                                            limit = request_limit,
                                            offset = request_offset,
                                            modified_before = access_date,
                                            format = "json"))

      } else {
        # filter gene set by title
        tribe_req <- httr::GET("https://tribe.greenelab.com/api/v1/geneset/",
                               query = list(title__startswith = type,
                                            show_versions = "true",
                                            organism = "9",  # P.a. is number 9
                                            xrdb = "Symbol",
                                            limit = request_limit,
                                            offset = request_offset,
                                            modified_before = access_date,
                                            format = "json"))
      }

      tribe_content <- suppressMessages(httr::content(tribe_req, as = "text"))
      gene_sets <- jsonlite::fromJSON(tribe_content)

      if (gene_sets$meta$total_count > 0) {
        gene_sets_df <- dplyr::data_frame(title = gene_sets$objects$title,
                                          genes = sapply(gene_sets$objects$versions,
                                                         function(x) unlist(x[1, "genes"])))
        gene_sets_df$count <- sapply(gene_sets_df$genes, function(x) length(x))
        gene_sets_df_list[[request_time]] <- gene_sets_df

        # the actual number of returned gene sets
        returned_size <- nrow(gene_sets_df)
        # If the return_size reaches the request limit, it indicates that we haven't
        # downloaded all gene sets yet. We will make a new request repeatly until
        # the number of returned gene sets is smaller than the request limit.
        if (returned_size < request_limit) {
          break
        } else {
          request_time <- request_time + 1
        }
      } else {
        break
      }
    }
  }

  if (length(gene_sets_df_list) == 0) {
    stop("Failed to retrieve any gene set. Please check the input parameters.")
  }

  # combine gene sets downloaded in all requests
  gene_sets_df <- dplyr::bind_rows(gene_sets_df_list)

  # subset gene sets using max_size and min_size
  gene_sets_df <- gene_sets_df[gene_sets_df$count >= min_size &
                                 gene_sets_df$count <= max_size, ]

  # extract gene sets list and name each gene set
  gene_sets_list <- gene_sets_df$genes
  names(gene_sets_list) <- gene_sets_df$title

  return(gene_sets_list)
}


#' Signatures annotation with gene sets
#'
#' Annotate signatures with known GO or KEGG gene sets downloaded from the Tribe
#' webserver. Annotation is done through a gene set enrichment test (fisher
#' exact test) with FDR correcting the number of gene sets tested for each
#' signature.
#'
#' @param selected_signatures a character vector storing names of signatures
#' @param model an ADAGE model to extract signatures from
#' @param genesets a named list of gene sets returned by the
#' function fetch_geneset(), each element in the list is a character vector
#' storing genes (PAO1 locus tag) in a gene set.
#' @param significance_cutoff numeric, FDR significance cutoff used to filter
#' gene sets (default: 0.05).
#' @return a data.frame storing significantly enriched gene sets for the
#' input signatures.
#' @export
annotate_signatures_with_genesets <- function(selected_signatures,
                                              model, genesets,
                                              significance_cutoff = 0.05){
  if (!check_input(model)) {
    stop("The model should be a data.frame with the first column as a character
         of gene IDs and the rest of the columns storing numeric weight values
         for each node.")
  }

  # extract all signatures' gene lists
  signatures_genes <- extract_signatures(model)
  # get genes in all signatures
  all_genes <- unique(unlist(signatures_genes))
  # make sure the input selected_signatures can be found in the model
  if (!all(selected_signatures %in% names(signatures_genes))){
    stop("Selected signatures can not be found in the specified model.")
  }
  # get the gene lists of only the selected signatures
  selected_signatures_genes <- signatures_genes[selected_signatures]

  # enrichment test between every signature - term pair
  signature_geneset <- lapply(seq_along(selected_signatures_genes), function(y){
    pvalue <- sapply(genesets, function(x)
      enrich_test(x, selected_signatures_genes[[y]], all_genes)$p.value)
    # correct by the number of gene sets tested
    qvalue <- p.adjust(pvalue, method = "fdr")
    # build a result table
    result_tb <- data.frame(signature = names(selected_signatures_genes)[y],
                            geneset = names(pvalue), pvalue = pvalue,
                            qvalue = qvalue, row.names = NULL)
    # filter by significance cutoff
    result_tb[result_tb$qvalue <= significance_cutoff, ]
  })
  signature_geneset <- dplyr::bind_rows(signature_geneset)

  return(signature_geneset)
}


#' Associated signatures retrieval
#'
#' Returns the significantly associated signatures for the input genesets. If
#' signature_limma_result is provided, it will only return the significantly
#' associated signatures that achieves the lowest adjusted p-value in the
#' signature_limma_result.
#'
#' @param input_genesets a character vector storing gene set names, must be
#' found in the provided gene set list
#' @param gene_set_list a named list storing all possible gene sets with each
#' element being a vector of gene IDs.
#' @param model an ADAGE model to extract signatures from
#' @param significance_cutoff numeric, FDR significance cutoff used to determine
#' signature-geneset association (default: 0.05).
#' @param signature_limma_result a data.frame that stores the result table
#' returned by limma. It includes logFC, adj.P.Val, and other statistics for
#' each signature.
#' @return a list storing associated signatures for each input gene set
#' @export
find_associated_signatures <- function(input_genesets, gene_set_list,
                                       model, significance_cutoff = 0.05,
                                       signature_limma_result = NULL){

  if(!all(input_genesets %in% names(gene_set_list))) {
    stop("The input gene sets cannot be found in the gene set list.")
  }

  # get the names of all signatures
  all_sigs <- c(paste0(colnames(model)[-1], "pos"),
                paste0(colnames(model)[-1], "neg"))

  # peform pathway association for all signatures
  all_enrichment <- annotate_signatures_with_genesets(
    selected_signatures = all_sigs, model = model, genesets = gene_set_list,
    significance_cutoff = significance_cutoff)

  # for each pathway, get its associated signatures
  associated_sigs <- lapply(input_genesets, function(x){
    sigs <- all_enrichment$signature[grepl(x, all_enrichment$geneset)]
    if (identical(sigs, character(0))) {
      # return NA is find no signature
      NA
    } else {
      if (is.null(signature_limma_result)) {
        sigs
      } else {
        # only record the one that is most significantly differentially active
        sigs[which.min(signature_limma_result[sigs, "adj.P.Val"])]
      }
    }
  })

  return(associated_sigs)
}

#' Annotating genes in signatures
#'
#' Annotates genes in the input signatures with their symbols, descriptions,
#' operons, associated signatures, and curated pathways if provided.
#'
#' @param selected_signatures a character vector storing names of signatures
#' @param model an ADAGE model to extract signatures from
#' @param curated_pathways a named list with each element being a gene set, such
#' as the output of the function fetch_geneset(). (default: NULL).
#' @return a data.frame storing genes in the input signatures. Each gene is
#' annotated by gene symbol, gene description, its operon, and signatures it
#' is in.
#' @export
annotate_genes_in_signatures <- function(selected_signatures, model,
                                         curated_pathways = NULL){

  if (!check_input(model)) {
    stop("The model should be a data.frame with the first column as a character
         of gene IDs and the rest of the columns storing numeric weight values
         for each node.")
  }

  # extract all signatures' gene lists
  signatures_genes <- extract_signatures(model)

  # make sure the input selected_signatures can be found in the model
  if (!all(selected_signatures %in% names(signatures_genes))){
    stop("Selected signatures can not be found in the specified model.")
  }

  # get the gene lists of only the selected signatures
  selected_signatures_genes <- signatures_genes[selected_signatures]

  # build a gene-signature map from the selected signatures
  gene_signature_map <- build_gene_geneset_map(selected_signatures_genes)
  colnames(gene_signature_map)[2] <- "signature"
  # replace ; with , because visNetwork can only use , to separate multiple groups
  gene_signature_map$signature <- gsub(";", ",", gene_signature_map$signature)
  gene_signature_map$N_signatures <- sapply(gene_signature_map$signature,
    function(x) length(unlist(strsplit(x,","))))

  # build a gene-operon map
  gene_operon_map <- build_gene_geneset_map(operons)
  colnames(gene_operon_map)[2] <- "operon"

  # incorporate gene symbol, description, operon, and signature annotations
  genes_df <- suppressWarnings(dplyr::right_join(gene_operon_map,
                                                 gene_signature_map))

  # annotate pathways each gene participate if curated_pathways is provided
  if (!is.null(curated_pathways)) {
    gene_pathway_map <- build_gene_geneset_map(curated_pathways)
    colnames(gene_pathway_map)[2] <- "pathway"
    genes_df <- suppressWarnings(dplyr::left_join(genes_df, gene_pathway_map))
  }

  genes_df <- suppressWarnings(
    dplyr::right_join(geneinfo[, c("LocusTag", "Symbol", "description")],
                                genes_df, by = c("LocusTag" = "geneID")))

  return(genes_df)

}


#' Gene-geneset map
#'
#' Map genes to the genesets they are in.
#'
#' @param genesets a named list with each element storing genes in one
#' geneset.
#' @return a data.frame with the first column specifying geneID and the second
#' column specifying the names of genesets that a gene is in.
build_gene_geneset_map <- function(genesets){

  # initialize a gene-geneset data.frame with a geneID column containing all
  # unique genes in the input genesets and a geneset column set to NA
  gene_geneset_df <- data.frame(geneID = unique(unlist(genesets)),
                                geneset = NA)

  # loop through each gene in each geneset
  for (set in names(genesets)) {
    for (gene in genesets[[set]]) {
      if (is.na(gene_geneset_df[gene_geneset_df$geneID == gene, "geneset"])){
        # if NA, this is the first gene-geneset relationship found for this
        # gene, replace NA with the geneset name
        gene_geneset_df[gene_geneset_df$geneID == gene, "geneset"] <- set
      }
      else{
        # this gene is already in some other genesets, simply append this
        # geneset's name to existing geneset names
        gene_geneset_df[gene_geneset_df$geneID == gene, "geneset"] <-
          paste(gene_geneset_df[gene_geneset_df$geneID == gene, "geneset"],
                set, sep = ";")
      }
    }
  }

  return(gene_geneset_df)
}
