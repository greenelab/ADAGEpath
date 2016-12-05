#' @importFrom magrittr "%>%"
#' @export
visualize_gene_networks <- function(selected_signatures, model = eADAGEmodel,
                                    cor_cutoff = 0.5, gene_logFC = NULL) {

  signatures_genes <- extract_signatures(model)
  selected_signatures_genes <- signatures_genes[selected_signatures]
  selected_genes <- unique(unlist(selected_signatures_genes))
  model <- model[model$geneID %in% selected_genes, ]

  # extract the weight matrix from the model
  weight_matrix <- as.matrix(model[, -1])
  rownames(weight_matrix) <- as.data.frame(model)[, 1]

  # calculate weight correlation
  weight_cor <- cor(t(weight_matrix))

  # everything below the cutoff is set to zero and will be ignored when
  # building the gene-gene network
  weight_cor[weight_cor < cor_cutoff & weight_cor > -cor_cutoff] <- 0

  # build an igraph network structure from the weight correlation matrix
  graph  <- igraph::graph_from_adjacency_matrix(weight_cor, weighted = TRUE,
                                                mode = 'undirected',
                                                diag = FALSE)

  # convert igraph structure to visNetwork structure
  gene_network <- visNetwork::toVisNetworkData(graph)

  # refine visualization options
  if (!is.null(gene_logFC)) {
    if (!is.null(selected_signatures)) {
      gene_logFC <- gene_logFC[
        gene_logFC$geneID %in% selected_genes, ]
    }
    gene_logFC$color <- create_logFC_colors(gene_logFC$logFC)
    gene_network$nodes <- suppressWarnings(dplyr::left_join(gene_network$nodes,
                                                            gene_logFC,
                                                            by = c("id" = "geneID")))
  }

  gene_network$nodes$label <- sapply(gene_network$nodes$label,
                                     function(x) suppressWarnings(to_symbol(x)))
  gene_network$nodes$size <- 15
  gene_signature_df <- gene_signature_map(selected_signatures_genes)
  gene_network$nodes <- suppressWarnings(dplyr::left_join(gene_network$nodes,
                                                          gene_signature_df,
                                                          by = c("id" = "gene")))
  gene_network$nodes$title <- paste0("<p>locus tag:", gene_network$nodes$id,
                                     "<br>symbol:", gene_network$nodes$label,
                                     "<br>logFC:", round(gene_network$nodes$logFC, 2),
                                     "<br>signatures:", gene_network$nodes$signature,
                                     "</p>")

  gene_network$edges$color <- ifelse(gene_network$edges$weight > 0,
                                     "#d01c8b", "#4dac26")
  min_weight <- min(abs(abs(gene_network$edges$weight)))
  range_weight <- diff(range(abs(gene_network$edges$weight)))
  gene_network$edges$width <- (abs(gene_network$edges$weight) - min_weight) /
    range_weight * 5
  gene_network$edges$shared_signatures <- two_genes_signature_map(
    gene_network$edges, selected_signatures_genes)
  gene_network$edges$title <- paste0("<p>correlation:",
                                     round(gene_network$edges$weight, 2),
                                     "<br>shared signatures:",
                                     gene_network$edges$shared_signatures,
                                     "</p>")

  final_network <- visNetwork::visNetwork(nodes = gene_network$nodes,
                   edges = gene_network$edges) %>%
        visNetwork::visPhysics(stabilization = FALSE) %>%
        visNetwork::visEdges(smooth = FALSE) %>%
        visNetwork::visOptions(selectedBy = list(variable = "signature",
                                                 multiple = TRUE),
                               highlightNearest = TRUE,
                               nodesIdSelection = TRUE) %>%
        visNetwork::visInteraction(navigationButtons = TRUE) %>%
        visNetwork::visLayout(randomSeed = 123)

  return(final_network)
}


create_logFC_colors <- function(logFC){
  col <- logFC
  pos_pal <- leaflet::colorNumeric(palette = "Reds",
                                   domain = logFC[logFC > 0])
  col[logFC > 0] <- pos_pal(logFC[logFC > 0])
  neg_pal <- leaflet::colorNumeric(palette = "Blues",
                                   domain = abs(logFC[logFC <= 0]))
  col[logFC <= 0] <- neg_pal(abs(logFC[logFC <= 0]))

  return(col)
}


gene_signature_map <- function(signatures_genes){
  gene_signature_df <- data.frame(gene = unique(unlist(signatures_genes)),
                                  signature = NA)
  for (sig in names(signatures_genes)) {
    for (gene in signatures_genes[[sig]]) {
      if (is.na(gene_signature_df[gene_signature_df$gene == gene, "signature"])){
        gene_signature_df[gene_signature_df$gene == gene, "signature"] <- sig
      }
      else{
        gene_signature_df[gene_signature_df$gene == gene, "signature"] <-
          paste(gene_signature_df[gene_signature_df$gene == gene, "signature"],
                sig, sep = ",")
      }
    }
  }

  return(gene_signature_df)

}

two_genes_signature_map <- function(gene_network_edge, signatures_genes){
  two_gene_sigs <- rep(NA, nrow(gene_network_edge))
  for (i in 1:nrow(gene_network_edge)) {
    have_sig_flag <- sapply(signatures_genes, function(x)
      all(gene_network_edge[i, c("from", "to")] %in% x))
    if (any(have_sig_flag)) {
      two_gene_sigs[i] <- paste(names(signatures_genes)[have_sig_flag], collapse = ",")
    }
  }

  return(two_gene_sigs)
}
