library("readr")

probedistribution_file <- "data-raw/compendium_quantile_distribution.txt"
probedistribution <- readr::read_lines(probedistribution_file)
probedistribution <- as.numeric(probedistribution)

devtools::use_data(probedistribution, overwrite = TRUE)
