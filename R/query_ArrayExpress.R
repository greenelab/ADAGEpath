
get_sample_annotation <- function(accession){

  if (!check_accession(accession)) {
    stop()
  } else {
    sample_info <- readr::read_tsv(
      paste0("https://www.ebi.ac.uk/arrayexpress/files/",
             accession, "/", accession, ".sdrf.txt"))
  }

  # move the "Array Data File" column to front
  sample_info <- sample_info[, c("Array Data File",
                                 setdiff(colnames(sample_info),
                                         "Array Data File"))]
  return(sample_info)

}


check_accession <- function(accession){
  AE_req <- httr::GET("https://www.ebi.ac.uk/arrayexpress/json/v3/experiments",
                      query = list(accession = accession))
  AE_content <- suppressWarnings(httr::content(AE_req, as = "text"))
  AE_exp <- jsonlite::fromJSON(AE_content)

  if (AE_exp$experiments$total == 0) {
    print("This experiment's accession number does not exist in ArrayExpress!")
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


download_raw <- function(accession, download_folder = "./"){

  AE_req <- httr::GET("https://www.ebi.ac.uk/arrayexpress/json/v3/files",
                      query = list(accession = accession))
  AE_content <- suppressWarnings(httr::content(AE_req, as = "text"))
  AE_file <- jsonlite::fromJSON(AE_content)
  all_files <- AE_file$files$experiment$file

  # check whether the raw file exists
  if (any(grepl("raw", all_files$name))) {

    # get the url of the raw data when it exists
    raw_url <- all_files$url[which(grepl("raw", all_files$name))]
    # specify the file path to save the downloaded the zip file
    raw_zip <- file.path(download_folder, basename(raw_url))
    # download it from the url
    download.file(raw_url, destfile = raw_zip)

    return(raw_zip)
  } else {

    stop("This experiment does not have raw cel files on ArrayExpress.")
  }
}
