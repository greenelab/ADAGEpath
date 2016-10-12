library("readr")

compendium_file <- "data-raw/all-pseudomonas-gene.pcl"
compendium <- read_tsv(compendium_file)
colnames(compendium)[1] <- "geneID"

devtools::use_data(compendium, overwrite = TRUE, compress = "bzip2")
