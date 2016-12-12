library("readr")
library("dplyr")

experiment_list <- readr::read_tsv("./data-raw/all-pseudomonas_SampleList.txt",
                                   col_names = FALSE)
experiment_list <- experiment_list[-1, ]

experimentID <- lapply(1:nrow(experiment_list), function(x){
  cbind(experiment_list[x, 1],
        unlist(strsplit(as.character(experiment_list[x, 2]), ";")))
})

experimentID <- dplyr::bind_rows(experimentID)
colnames(experimentID) <- c("Experiment", "Sample")

devtools::use_data(experimentID, overwrite = TRUE, compress = "bzip2")
