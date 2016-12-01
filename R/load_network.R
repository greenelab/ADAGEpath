#' @importFrom magrittr "%>%"
visualize_gene_networks <- function(model = eADAGEmodel, cor_cutoff = 0.5,
                                    selected_signatures = NULL,
                                    gene_logFC = NULL){

  if (!is.null(selected_signatures)){
    all_signatures <- extract_signatures(model)
    selected_genes <- unique(unlist(all_signatures[selected_signatures]))
    model <- model[model$geneID %in% selected_genes, ]
  }

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

  if (!is.null(gene_logFC)) {
    if (!is.null(selected_signatures)) {
      gene_logFC <- gene_logFC[
        gene_logFC$geneID %in% selected_genes, ]
    }
    gene_logFC$color <- create_logFC_colors(gene_logFC$logFC)
    gene_network$nodes$color <- gene_logFC$color[
      match(gene_network$nodes$label, gene_logFC$geneID)]
  }

  gene_network$nodes$label <- sapply(gene_network$nodes$label,
                                     function(x) to_symbol(x))

  gene_network$edges$width <- abs(gene_network$edges$weight*3)
  gene_network$edges$color <- ifelse(gene_network$edges$weight > 0,
                                     "blue", "yellow")

  print(visNetwork::visNetwork(nodes = gene_network$nodes,
                   edges = gene_network$edges) %>%
        visNetwork::visPhysics(stabilization = FALSE) %>%
        visNetwork::visEdges(smooth = FALSE))

}


create_logFC_colors <- function(logFC){
  col <- logFC
  pos_pal <- leaflet::colorNumeric(palette = "Reds",
                                   domain = logFC[logFC > 0])
  col[logFC > 0] <- pos_pal(logFC[logFC > 0])
  neg_pal <- leaflet::colorNumeric(palette = "Greens",
                                   domain = abs(logFC[logFC <= 0]))
  col[logFC <= 0] <- neg_pal(abs(logFC[logFC <= 0]))

  return(col)
}


