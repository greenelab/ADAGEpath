library("readr")

compendium_file <- "data-raw/all-pseudomonas-gene.pcl"
PAcompendium <- readr::read_tsv(compendium_file)
colnames(PAcompendium)[1] <- "geneID"

devtools::use_data(PAcompendium, overwrite = TRUE, compress = "bzip2")
