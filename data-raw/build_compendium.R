library("httr")
library("jsonlite")
library("affy")
library("affyio")
library("Biobase")
library("preprocessCore")
library("tibble")

# This script prepares compendium-related data objects for the package. They
# include "PAcompendium", "probedistribution", and "experimentID". It
# downloads datasets from the ArrayExpress database, processes the raw cel
# files into an expression compendium, builds a map between experiment
# accession number and cel file, and records the quantile distribution at the
# probe level. The current setting allows it to build the same compendium as
# used in the eADAGE paper. Modify "end_date" and "exclude_samples" to build an
# up-to-date compendium.

# set parameters--------------

# parameters used to query ArrayExpress
species <- "Pseudomonas aeruginosa"
array <- "A-AFFY-30"
raw <- "true"
start_date <- "1990-01-01"
end_date <- "2015-07-31"
# only download datasets released in this date range
date <- paste0("[", start_date, " ", end_date, "]")

# To precisely re-build the compendium used in the eADAGE paper, the following
# two samples are excluded because they were added after 2015-07-31, though
# their experiment was released before 2015-07-31.
exclude_samples <- c("aerobic_NO3_1.CEL", "aerobic_NO3_2.CEL")

# folders to save downloaded files
download_folder <- "./data-raw/zips"
cel_folder <- "./data-raw/cels"
dir.create(download_folder)
dir.create(cel_folder)

# file path to save the built compendium
compendium_file <- paste0("./data-raw/PAcompendium_", end_date, ".txt")
chip_type <- "Pae_G1a"

# query ArrayExpress--------------

AE_req <- httr::GET("https://www.ebi.ac.uk/arrayexpress/json/v3/files",
                    query = list(species = species,
                                 array = array,
                                 raw = raw,
                                 date = date))
AE_content <- suppressWarnings(httr::content(AE_req, as = "text"))
AE_files <- jsonlite::fromJSON(AE_content)
N_experiments <- AE_files$files$`total-experiments`

# download datasets--------------

experiment_sample <- vector("list", N_experiments)
for (i in 1:N_experiments) {
  accession <- AE_files$files$experiment$accession[[i]]
  experiment_files <- AE_files$files$experiment$file[[i]]
  raw_url <- experiment_files[which(experiment_files$kind == "raw"), "url"]
  raw_zip <- file.path(download_folder, basename(raw_url))
  download.file(raw_url, destfile = raw_zip)
  cel_files <- unzip(raw_zip, exdir = cel_folder)
  experiment_sample[[i]] <- data.frame(Experiment = accession,
                                       Sample = basename(cel_files),
                                       stringsAsFactors = FALSE)
}

# process cel files into expression values---------------

# list all celfiles in the directory to be processed
celfiles <- affy::list.celfiles(cel_folder)
# exclude some celfiles
celfiles <- setdiff(celfiles, exclude_samples)
# ptype is now a vector with the type of each array
ptype <- sapply(celfiles,
                function(f) affyio::read.celfile.header(paste(cel_folder, f,
                                                              sep = "/"))[1])
# pfiles vector only contains celfiles of the specified chip type
pfiles <-  subset(celfiles, ptype == chip_type)
pfiles_path <- file.path(cel_folder, pfiles)
# ReadAffy loads the array data using the custom CDF
affy_batch <- affy::ReadAffy(filenames = pfiles_path)
# rma processes the data by multi-array average expression measure
expressionSet <- affy::rma(affy_batch)
# convert expressionSet to a data.frame
PAcompendium <- data.frame(Biobase::exprs(expressionSet),
                           stringsAsFactors = FALSE, check.names = FALSE)
PAcompendium <- PAcompendium[, order(colnames(PAcompendium))]
PAcompendium <- tibble::rownames_to_column(PAcompendium, var = "geneID")
# only preserve gene IDs that are PA numbers
PAcompendium <- PAcompendium[startsWith(PAcompendium$geneID, "PA"), ]
PAcompendium$geneID <- sapply(PAcompendium$geneID,
                              function(x) unlist(strsplit(x, "_"))[1])
rownames(PAcompendium) <- NULL
# save the compendium
write.table(PAcompendium, compendium_file, sep = "\t", quote = FALSE,
            row.names = FALSE)
devtools::use_data(PAcompendium, overwrite = TRUE, compress = "bzip2")

# build an experiment-sample map-----------------

experimentID <- do.call(rbind, experiment_sample)
experimentID <- experimentID[experimentID$Sample %in% pfiles, ]
devtools::use_data(experimentID, overwrite = TRUE, compress = "bzip2")

# record the quantile distribution at the probe level---------------

# return perfect match probes
PMmat <- affy::pm(affy_batch, NULL)
# background correction
PMmat_bg <- preprocessCore::rma.background.correct(PMmat)
# calculate the quantile distribution
probedistribution <- preprocessCore::normalize.quantiles.determine.target(
  PMmat_bg)
devtools::use_data(probedistribution, overwrite = TRUE)
