#' Gene set fetch
#'
#' Fetches GO or KEGG gene sets from the Tribe webserver
#' (url{tribe.greenelab.com}).
#'
#' @param type character, type of the gene sets, either "GO" or "KEGG"
#' (default: "KEGG").
#' @param max_size int, maximum gene set size to be considered as a meaningful
#' gene set (default: 100).
#' @param min_size int, minimum gene set size to be considered as a meaningful
#' gene set (default: 5).
#' @return a named list with each element being a gene set
#' @export
fetch_geneset <- function(type = "KEGG", max_size = 100, min_size = 5){

  # make sure type is one of "GO" or "KEGG"
  if (!type %in% c("GO", "KEGG")) {
    stop("type can only be either GO or KEGG.")
  }

  # to prevent overloading the TRIBE webserver, we limit the download size
  # to 2000 gene sets at one time and make multiple downloads through changing
  # the offset until all gene sets have been downloaded.
  request_limit <- 2000
  request_time <- 1
  gene_sets_df_list <- list()
  while (TRUE) {
    # set request offset to the number of request limits that have already
    # been reached
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
    # If the return_size reaches the request limit, it indicates that we haven't
    # downloaded all gene sets yet. We will make a new request repeatly until
    # the number of returned gene sets is smaller than the request limit.
    if (returned_size < request_limit) {
      break
    } else {
      request_time <- request_time + 1
    }
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
#' (default: the 300-node eADAGE model preloaded in the package).
#' @param genesets a named list of gene sets returned by the
#' function fetch_geneset(), each element in the list is a character vector
#' storing genes (PAO1 locus tag) in a gene set.
#' @param significance_cutoff numeric, FDR significance cutoff used to filter
#' the result (default: 0.05).
#' @return a data.frame storing significantly enriched gene sets for the
#' input signatures.
#' @export
annotate_signatures_with_genesets <- function(selected_signatures,
                                              model = eADAGEmodel, genesets,
                                              significance_cutoff = 0.05){
  if (!check_input(model)) {
    stop("The model should be a data.frame with first column being gene IDs
         in character and the rest columns storing numeric weight values for
         each node per column.")
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


#' Annotating genes in signatures
#'
#' Annotates genes in the input signatures with their symbols, descriptions,
#' associated signatures.
#'
#' @param selected_signatures a character vector storing names of signatures
#' @param model an ADAGE model to extract signatures from
#' (default: the 300-node eADAGE model preloaded in the package).
#' @return a data.frame storing genes in the input signatures. Each gene is
#' annotated by gene symbol, gene description, and signatures it is in.
#' @export
annotate_genes_in_signatures <- function(selected_signatures,
                                         model = eADAGEmodel){

  if (!check_input(model)) {
    stop("The model should be a data.frame with first column being gene IDs
         in character and the rest columns storing numeric weight values for
         each node per column.")
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
  gene_signature_map <- build_gene_signature_map(selected_signatures_genes)
  gene_signature_map$N_signatures <- sapply(gene_signature_map$signature,
    function(x) length(unlist(strsplit(x,","))))

  # build a gene-operon map
  gene_operon_map <- build_gene_signature_map(operons)
  colnames(gene_operon_map)[2] <- "operon"

  # incorporate gene symbol, description, operon, and signature annotations
  genes_df <- suppressWarnings(dplyr::right_join(gene_operon_map,
                                                 gene_signature_map))
  genes_df <- suppressWarnings(
    dplyr::right_join(geneinfo[, c("LocusTag", "Symbol", "description")],
                                genes_df, by = c("LocusTag" = "geneID")))

  return(genes_df)

}


#' Gene-signature map
#'
#' Map genes to the signatures they are in.
#'
#' @param signatures_genes a named list with each element storing genes in one
#' signature.
#' @return a data.frame with the first column specifying geneID and the second
#' column specifying signatures that a gene is in.
build_gene_signature_map <- function(signatures_genes){

  # initialize a gene-signature data.frame with a geneID column containing all
  # unique genes in the input signatures and a signature column set to NA
  gene_signature_df <- data.frame(geneID = unique(unlist(signatures_genes)),
                                  signature = NA)

  # loop through each gene in each signature
  for (sig in names(signatures_genes)) {
    for (gene in signatures_genes[[sig]]) {
      if (is.na(gene_signature_df[gene_signature_df$geneID == gene, "signature"])){
        # if NA, this is the first gene-signature relationship found for this
        # gene, replace NA with the signature name
        gene_signature_df[gene_signature_df$geneID == gene, "signature"] <- sig
      }
      else{
        # this gene is already in some other signatures, simply append this
        # signature's name to existing signature names
        gene_signature_df[gene_signature_df$geneID == gene, "signature"] <-
          paste(gene_signature_df[gene_signature_df$geneID == gene, "signature"],
                sig, sep = ",")
      }
    }
  }

  return(gene_signature_df)
}
