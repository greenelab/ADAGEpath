#' Gene-gene network visualization
#'
#' Builds and visualizes an ADAGE-derived gene-gene network of genes in the
#' selected signatures. In the network, each node represents a gene and two
#' nodes are linked by an edge if the correlation of the two connected genes in
#' the ADAGE model is higher than a certain cutoff. Edge width denotes the
#' correlation strength and edge color indicates positive (magenta) or negative
#' (green) correlation. If a value associated with the expression of each gene in
#' a dataset (such as fold change or correlation to phenotype) is provided, it
#' will be used to color gene nodes in the network with red meaning high positive
#' value and blue meaning high negative value.
#'
#' @param selected_signatures a vector storing names of signatures whose genes
#' to include in the gene-gene network
#' @param model an ADAGE model to build gene-gene network from
#' (default: the 300-node eADAGE model preloaded in the package).
#' @param cor_cutoff numeric, the correlation cutoff to decide whether an edge
#' between two genes exists in the network (default: 0.5).
#' @param gene_color_value a data.frame with two columns, one is geneID,
#' the other could be logFC, correlation, or other value associated with
#' each gene. This value will be used to color genes in the network.
#' If not provided, the gene color in the gene-gene network is uniformly blue.
#' (default:NULL).
#' @importFrom magrittr "%>%"
#' @export
visualize_gene_network <- function(selected_signatures, model = eADAGEmodel,
                                    cor_cutoff = 0.5, gene_color_value = NULL) {

  # extract all signatures and their genes from the model using the custom
  # function extract_signatures()
  signatures_genes <- extract_signatures(model)

  # make sure the input selected_signatures can be found in the model
  if (!all(selected_signatures %in% names(signatures_genes))){
    stop("Names of the selected signatures are not in the specified ADAGE model.")
  }

  # only include genes in the selected signatures
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

  # annotate genes in the selected signatures
  gene_annotation <- annotate_genes_in_signatures(selected_signatures)
  gene_network$nodes <- suppressWarnings(dplyr::left_join(gene_network$nodes,
                                                          gene_annotation,
                                                          by = c("id" = "LocusTag")))
  gene_network$nodes <- gene_network$nodes %>% dplyr::select(-label) %>%
    dplyr::rename(label = Symbol)

  # set network node size
  gene_network$nodes$size <- 15

  if (!is.null(gene_color_value)) {
    # set network node color to reflect gene fold change
    gene_color_value <- gene_color_value[gene_color_value$geneID %in%
                                         selected_genes, ]
    value_name <- colnames(gene_color_value)[2]
    gene_color_value$color <- map_value_to_color(gene_color_value[, 2])
    gene_network$nodes <- suppressWarnings(dplyr::left_join(gene_network$nodes,
                                                            gene_color_value,
                                                            by = c("id" = "geneID")))

    # set network node title that will be displayed when mouse is above the node
    gene_network$nodes$title <- paste0("<p>locus tag:", gene_network$nodes$id,
                                       "<br>symbol:", gene_network$nodes$label,
                                       "<br>operon:", gene_network$nodes$operon,
                                       "<br>description:", gene_network$nodes$description,
                                       "<br>",value_name, ":",
                                       round(gene_network$nodes[, value_name], 2),
                                       "<br>signatures:", gene_network$nodes$signature,
                                       "</p>")
  } else {
    # set network node title that will be displayed when mouse is above the node
    # (without fold change)
    gene_network$nodes$title <- paste0("<p>locus tag:", gene_network$nodes$id,
                                       "<br>symbol:", gene_network$nodes$label,
                                       "<br>operon:", gene_network$nodes$operon,
                                       "<br>description:", gene_network$nodes$description,
                                       "<br>signatures:", gene_network$nodes$signature,
                                       "</p>")
  }

  # set network edge color to magenta if edge weight (correlation) is positive
  # and green if negative
  gene_network$edges$color <- ifelse(gene_network$edges$weight > 0,
                                     "#d01c8b", "#4dac26")

  # set the network edge width to be linear to edge weight
  min_weight <- min(abs(abs(gene_network$edges$weight)))
  range_weight <- diff(range(abs(gene_network$edges$weight)))
  gene_network$edges$width <- (abs(gene_network$edges$weight) - min_weight) /
    range_weight * 5

  # annotate network edge to include a column indicating which signatures two
  # genes share with each other
  gene_network$edges <- annotate_shared_signatures(
    gene_network$edges, selected_signatures_genes)

  # set network edge title that will be displayed when mouse is above the edge
  gene_network$edges$title <- paste0("<p>correlation:",
                                     round(gene_network$edges$weight, 2),
                                     "<br>shared signatures:",
                                     gene_network$edges$shared_signatures,
                                     "</p>")

  # set the final network display options
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


#' Mapping value to HEX color
#'
#' Maps the continuous values to HEX color values. The positive values
#' are mapped to the "Reds" palette from Color Brewer 2
#' \url{http://colorbrewer2.org/#type=sequential&scheme=Reds&n=3},
#' and the negative values
#' are oppsitely mapped to the "Blues" palette from the Color Brewer 2
#' \url{http://colorbrewer2.org/#type=sequential&scheme=Blues&n=3}.
#'
#' @param value a numeric vector storing values to map
#' @return a character vector storing the mapped HEX color values
map_value_to_color <- function(value){
  col <- value
  pos_pal <- leaflet::colorNumeric(palette = "Reds",
                                   domain = value[value > 0])
  col[value > 0] <- pos_pal(value[value > 0])
  neg_pal <- leaflet::colorNumeric(palette = "Blues",
                                   domain = abs(value[value <= 0]))
  col[value <= 0] <- neg_pal(abs(value[value <= 0]))

  return(col)
}


#' Shared signature annotation
#'
#' Annotates shared signatures for gene pairs in the input gene_network_edge
#' data.frame.
#'
#' @param gene_network_edge a data.frame specifying gene-gene connections
#' in a gene-gene network. The first two columns are "to" and "from" geneIDs.
#' It is an intermediate product from the visualize_gene_network() function.
#' @param signatures_genes a named list with each element storing genes in one
#' signature.
#' @return the data.frame gene_network_edge with one more column
#' shared_signatures that specifies the shared signatures of the from and to
#' genes.
annotate_shared_signatures <- function(gene_network_edge, signatures_genes){

  # initialize the shared signature vector with NA
  two_gene_sigs <- rep(NA, nrow(gene_network_edge))

  # loop through each row of gene_network_edge that stores one from-to
  # gene pair
  for (i in 1:nrow(gene_network_edge)) {

    # find out all the signatures that both have "from" gene and "to" gene
    have_sig_flag <- sapply(signatures_genes, function(x)
      all(gene_network_edge[i, c("from", "to")] %in% x))

    if (any(have_sig_flag)) {
      # if there is one or more such signatures, paste their names and update
      # to return vector
      two_gene_sigs[i] <- paste(names(signatures_genes)[have_sig_flag],
                                collapse = ",")
    }
  }
  gene_network_edge$shared_signatures <- two_gene_sigs

  return(gene_network_edge)
}
