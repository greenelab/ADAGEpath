#' Querying sample information from ArrayExpress
#'
#' Gets sample information of an experiment stored in the SDRF (Sample and Data
#' Relationship Format) file from ArrayExpress database.
#'
#' @param accession ArrayExpress experiment accession number in the format of
#' "E-XXXX-n".
#' @return a data.frame storing sample information in the SDRF file
#' @export
get_sample_info <- function(accession){
  AE_req <- httr::GET("https://www.ebi.ac.uk/arrayexpress/json/v3/files",
                      query = list(accession = accession))
  AE_content <- suppressWarnings(httr::content(AE_req, as = "text"))
  AE_file <- jsonlite::fromJSON(AE_content)

  # check whether this accession number has experiments
  if (AE_file$files$`total-experiments` > 0) {
    all_files <- AE_file$files$experiment$file

    # check whether the sdrf file exists
    if (any(grepl("sdrf", all_files$name))) {
      # get the url of the sdrf when it exists
      sdrf_url <- all_files$url[which(grepl("sdrf", all_files$name))]
      # read in the sdrf file which is a tab-delimited file
      sample_info <- readr::read_tsv(sdrf_url)
      # move the "Array Data File" column to front
      sample_info <- sample_info[, c("Array Data File",
                                     setdiff(colnames(sample_info),
                                             "Array Data File"))]
      return(sample_info)
    } else {
      stop("This experiment does not have the sdrf file on ArrayExpress.")
    }
  } else {
    stop("This accession number has no associated experiments.")
  }
}


#' Accession number validity check
#'
#' Checks whether the accession number has an associated experiment and the
#' experiment was measured on the Pseudomonas affymetric platform "A-AFFY-30".
#'
#' @param accession ArrayExpress experiment accession number in the format of
#' "E-XXXX-n".
#' @return logical, TRUE if the accession number has an associated experiment
#' and the experiment is on the "A-AFFY-30" platform; otherwise FALSE.
check_accession <- function(accession){
  AE_req <- httr::GET("https://www.ebi.ac.uk/arrayexpress/json/v3/experiments",
                      query = list(accession = accession))
  AE_content <- suppressWarnings(httr::content(AE_req, as = "text"))
  AE_exp <- jsonlite::fromJSON(AE_content)

  if (AE_exp$experiments$total == 0) {
    print("This accession number has no associated experiment!")
    return(FALSE)
  } else if (AE_exp$experiments$experiment$arraydesign[[1]]$accession !=
             "A-AFFY-30") {
    print("This experiment is not performed on the P.a. Affy platform
          A-AFFY-30.")
    return(FALSE)
  } else {
    return(TRUE)
  }
}


#' Downloading raw data from ArrayExpress
#'
#' Downloads the raw cel files of an experiment from the ArrayExpress database.
#'
#' @param accession ArrayExpress experiment accession number in the format of
#' "E-XXXX-n".
#' @return file path to the downloaded zip file that stores the cel files of
#' the input experiment.
download_raw <- function(accession, download_folder = "./"){

  AE_req <- httr::GET("https://www.ebi.ac.uk/arrayexpress/json/v3/files",
                      query = list(accession = accession))
  AE_content <- suppressWarnings(httr::content(AE_req, as = "text"))
  AE_file <- jsonlite::fromJSON(AE_content)

  # check whether this accession number has experiments
  if (AE_file$files$`total-experiments` > 0) {
    all_files <- AE_file$files$experiment$file

    # check whether the raw file exists
    if (any(grepl("raw", all_files$name))) {

      # get the url of the raw data when it exists
      raw_url <- all_files$url[which(grepl("raw", all_files$name))]
      # create the download folder
      dir.create(download_folder)
      # specify the file path to save the downloaded the zip file
      raw_zip <- file.path(download_folder, basename(raw_url))
      # download it from the url
      download.file(raw_url, destfile = raw_zip)

      return(raw_zip)
    } else {
      stop("This experiment does not have raw cel files on ArrayExpress.")
    }
  } else {
    stop("This accession number has no associated experiments.")
  }
}
