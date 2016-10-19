library("readr")
library("dplyr")

compendium_file <- "data-raw/all-pseudomonas-gene.pcl"
compendium <- readr::read_tsv(compendium_file)
colnames(compendium)[1] <- "geneID"
geneID <- compendium[, 1]

networkfile <- "data-raw/net300_100models_1_100_k=300_seed=123_ClusterByweighted_avgweight_network_ADAGE.txt"
network <- readr::read_delim(networkfile, delim = "\t", col_names = F,
                             n_max = nrow(geneID), skip = 2)
colnames(network) <- paste0("Node", seq(1, ncol(network)))
eADAGEmodel <- dplyr::bind_cols(geneID, network)

devtools::use_data(eADAGEmodel, overwrite = TRUE, compress = "bzip2")
