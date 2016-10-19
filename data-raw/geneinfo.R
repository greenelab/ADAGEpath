library("readr")

ftp_url <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Archaea_Bacteria/"
gene_info_file <- "Pseudomonas_aeruginosa_PAO1.gene_info.gz"

download.file(paste0(ftp_url, gene_info_file),
              file.path("data-raw", gene_info_file))
geneinfo <- readr::read_tsv(file.path("data-raw", gene_info_file))

devtools::use_data(geneinfo, overwrite = TRUE, compress = "bzip2")
