library("readr")
library("dplyr")

compendium_file <- "data-raw/all-pseudomonas-gene.pcl"
compendium <- read_tsv(compendium_file)
colnames(compendium)[1] <- "geneID"

gene_max <- apply(compendium[, -1], 1, max)
gene_min <- apply(compendium[, -1], 1, min)
expressionrange <- data_frame(geneID = compendium[[1]],
                              max_express = gene_max,
                              min_express = gene_min)

devtools::use_data(expressionrange, overwrite = TRUE)
