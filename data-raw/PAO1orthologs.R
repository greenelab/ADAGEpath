library("readr")

file_url <- "http://pseudomonas.com/downloads/pseudomonas/pgd_r_16_1/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107_orthologs.txt"

PAO1orthologs <- readr::read_tsv(file_url)

devtools::use_data(PAO1orthologs, overwrite = TRUE, compress = "bzip2")
