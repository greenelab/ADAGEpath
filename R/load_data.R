#' Loading and processing input data
#'
#' Currently supports five types of inputs:
#' file path to a zip file containing microarray CEL files in a dataset
#' (isProcessed=FALSE, isRNAseq=FALSE);
#' file path to a folder containing microarray CEL files in a dataset
#' (isProcessed=FALSE, isRNAseq=FALSE);
#' an ArrayExpress experiment accession number in the format of "E-XXXX-n"
#' (the experiment must be measured on the "A-AFFY-30" platform)
#' (isProcessed=FALSE, isRNAseq=FALSE);
#' file path to a tab-delimited txt file storing processed gene expression
#' values in a dataset (gene identifiers in the first column and then one
#' sample per column) (isProcessed=TRUE, isRNAseq=TRUE/FALSE);
#' a data.frame object storing processed gene expression
#' values in a dataset (gene identifiers in the first column and then one
#' sample per column) (isProcessed=TRUE, isRNAseq=TRUE/FALSE);
#'
#' @param input file path to the input file or input folder or ArrayExpress
#' accession number or a data.frame object
#' @param isProcessed a logical value indicating whether the input_data has
#' already been processed into expression values at the gene level.
#' @param isRNAseq a logical value indicating whether the processed input_data
#' is RNAseq data. If TRUE, the processed RNAseq data will be normalized to
#' a comparable range with the microarray-based compendium using TDM. If FALSE,
#' the processed input_data is considered a microarray dataset  and will be
#' quantile normalized to be comparable to the compendium.
#' @param model the ADAGE model used to analyze the input_data
#' @param compendium the gene expression compendium of an organism
#' @param quantile_ref a vector storing the reference quantile distribution of
#' the input compendium at the microarray probe level.
#' Since the input microarray data needs to be normalized to the processed
#' compendium, the compendium and the quantile_ref must match each other.
#' @param download_folder file path to save files downloaded from ArrayExpress
#' when input is an ArrayExpress accession number.
#' @param norm01 a logical value indicating whether the output should be
#' zero-one normalized (default: FALSE)
#' @return a data.frame containing the processed gene expression values ready
#' for ADAGE analysis.
#' @export
load_dataset <- function(input, isProcessed, isRNAseq, model, compendium,
                         quantile_ref, download_folder = "./",
                         norm01 = FALSE){

  if (!check_input(model)) {
    stop("The model should be a data.frame with the first column as a character
         of gene IDs and the rest of the columns storing numeric weight values
         for each node.")
  }

  if (!check_input(compendium)) {
    stop("The compendium should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric gene
         expression values for each sample per column.")
  }

  # quantile normalize probes if the input is raw CEL format, directly load data
  # if the input is processed microarray or RNAseq data
  if (!isProcessed) {
    # if the input is a zip file, unzip into the same folder and use cel files
    # in the folder. If the input is a folder, directly use cel files in the
    # folder. If the input is ArrayExpress accession number, download its raw
    # data zip, upzip it, and use its cel files.

    if (endsWith(input, ".zip")) {
      # unzip files
      cel_folder <- file.path(dirname(input), gsub(".zip", "", basename(input)))
      unzip(input, exdir = cel_folder)


    } else if (R.utils::isDirectory(input)) {
      cel_folder <- input

    } else if (grepl("E-\\w{4}-\\d+", input)) {
      # check whether the accession number exists on ArrayExpress and the
      # associated experiment is measured on A-AFFY-30 platform
      if (!check_accession(input)) {
        stop()
      } else{
        raw_zip <- download_raw(input, download_folder = download_folder)
        cel_folder <- file.path(dirname(raw_zip), gsub(".zip", "",
                                                       basename(raw_zip)))
        unzip(raw_zip, exdir = cel_folder)

      }

    } else {
      stop("When isProcessed is FALSE, the input should be a zip file with CEL
           files inside, or a folder of CEL files, or a valid ArrayExpress
           accession number.")
    }

    # normalize each probe
    data <- process_celfiles(cel_folder = cel_folder, use_ref = TRUE,
                             quantile_ref = quantile_ref)

    # transform gene features in the input_data to gene features used in ADAGE model
    data <- match_IDs(input_data = data, ref_IDs = as.data.frame(model)[, 1])

  } else {

    if (is.data.frame(input)) {
      # expression data has been read into the data.frame input already
      data <- input
    } else if (file.exists(input)) {
      # read in the processed data when input is a file path
      data <- readr::read_tsv(input)
    } else {
      stop("When isProcessed is TRUE, the input should be an existing data.frame
           storing processed expression values or a file path to a tab-delimited
           file storing processed expression values.")
    }
    if (!check_input(data)) {
      stop("The input data should be a data.frame with the first column as a
           character of gene IDs and the rest of the columns storing numeric
           expression values for each sample.")
    }
    colnames(data)[1] <- "geneID"

    # transform gene features in the input_data to gene features used in ADAGE model
    data <- match_IDs(input_data = data, ref_IDs = as.data.frame(model)[, 1])

    # impute missing values if exist
    if (any(is.na(data))) {
      data <- impute_miss_values(input_data = data, ref_data = compendium)
    }

    if (isRNAseq) {

      # perform TDM transformation if the input is RNAseq data
      data <- TDM_RNAseq(input_data = data, ref_data = compendium)
    } else {

      # perform quantile normalization if the input is microarray data
      data <- quantile_norm(input_data = data, use_ref = TRUE,
                            ref_data = compendium)
    }
  }

  if (norm01) {
    # perform zero-one normalization
    data <- zeroone_norm(input_data = data, use_ref = TRUE,
                         ref_data = compendium)
  }

  # rownames are not needed as geneID is stored in the first column
  rownames(data) <- NULL

  return(data)
}


#' Processing CEL files
#'
#' Processes microarray data in CEL format. Only CEL files measured on the
#' "Pae_G1a" platform will be processed. If using reference, the quantile
#' normalization of probes is based on a reference quantile
#' distribution.
#'
#' @param cel_folder file path to the folder storing CEL files.
#' @param use_ref a logical value indicating whether the probe normalization
#' should be done using a reference quantile distribution (default: TRUE).
#' @param quantile_ref a vector storing the reference quantile distribution
#' @return a data.frame containing normalized gene expression values with geneID
#' in the first column and then each CEL file in one column.
process_celfiles <- function(cel_folder, use_ref = TRUE, quantile_ref){

  # read in all CEL files in the folder
  celfiles <- affy::list.celfiles(cel_folder)

  # get the platform types
  ptype <- sapply(celfiles, function(f) {
    affyio::read.celfile.header(file.path(cel_folder, f))[1]
  })

  # only process CEL files on "Pae_G1a" platform
  pfiles <- paste(cel_folder, subset(celfiles, ptype == "Pae_G1a"), sep = "/")

  # read in CEL files and produce an AffyBatch object
  affybatch <- affy::ReadAffy(filenames = pfiles)


  if (use_ref) {

    # return perfect match probes
    PMmat <- affy::pm(affybatch, NULL)

    # background correction
    PMmat_bg <- preprocessCore::rma.background.correct(PMmat)

    # quantile normalization using the reference distribution
    PMmat_normed <- preprocessCore::normalize.quantiles.use.target(
      PMmat_bg, target = quantile_ref)

    # get probe sets
    probe_list <- affy::probeNames(affybatch, NULL)

    # summarize probes into genes
    expression <- preprocessCore::subColSummarizeMedianpolishLog(
      PMmat_normed, group.labels = probe_list)

  } else {

    # directly use rma function to normalize within the input_data
    expression <- affy::rma(affybatch)

  }

  # build the final gene expression data.frame
  colnames(expression) <- celfiles
  expression <- data.frame(geneID = rownames(expression), expression,
                           stringsAsFactors = FALSE, check.names = FALSE)

  return(expression)

}


#' Converting gene ID to locus tag
#'
#' Converts the input gene ID to PAO1 locus tag that are used in ADAGE models.
#' It currently can convert gene IDs on AFFY chip, PAO1 gene symbols, and
#' PAO1 gene orthologs in other P.a. strains. It returns NA if the input is not
#' recognized as one of above.
#'
#' @param input_ID character, the input gene ID.
#' @param ref_IDs a vector storing reference gene IDs that do not need
#' conversion
#' @return the corresponding PAO1 locus tag ("PAXXXX") for the input gene or NA
#' if the input gene is not recognized.
#' @export
to_LocusTag <- function(input_ID, ref_IDs) {

  if (input_ID %in% ref_IDs) {

    # do nothing if the input ID is already in the ref_IDs
    return(input_ID)

  } else if (startsWith(input_ID, "ig") | startsWith(input_ID, "Pae") |
      startsWith(input_ID, "AFFY")) {

    # gene IDs start with "ig", "Pae", or "AFFY" are controls on the AFFY chip
    return("control")

  } else if (endsWith(input_ID, "_at")) {

    # gene IDs that end with "_at" but are not controls are AFFY IDs in the
    # format "PAXXXX_symbol", only the first part before "_" will be preserved
    return(unlist(strsplit(input_ID, "_"))[1])

  } else if (input_ID %in% geneinfo$Symbol) {

    # convert gene symbol to locus tag
    output_ID <- geneinfo$LocusTag[geneinfo$Symbol == input_ID]
    return(output_ID)

  } else if (input_ID %in% PAO1orthologs$`Locus Tag (Hit)`) {

    # map ortholog gene to PAO1 locus tag
    output_ID <- PAO1orthologs$`Locus Tag (Query)`[
      PAO1orthologs$`Locus Tag (Hit)` == input_ID]
    return(output_ID)

  } else {

    warning(paste("Gene", input_ID, "from the input file is not found in the
                  gene database and cannot be converted to PAO1 locus tag!"))
    return(NA)

  }
}


#' Converting locus tag to symbol
#'
#' Converts PAO1 gene locus tag to gene symbol. The input ID should be a
#' PAO1 locus tag. If it is not found in the PAO1 geneinfo database,
#' it will not be converted and the return value will be the same as the input.
#'
#' @param input_ID character, a PAO1 gene locus tag, such as "PA0001".
#' @return the converted gene symbol for the input locus tag.
#' @export
to_symbol <- function(input_ID){

  if (input_ID %in% geneinfo$Symbol) {
    return(input_ID)
  } else if (input_ID %in% geneinfo$LocusTag) {
    return(geneinfo$Symbol[geneinfo$LocusTag == input_ID])
  } else {
    warning(paste(input_ID, "not recognized as PAO1 locus tag.",
                  "It is kept unchanged."))
    return(input_ID)
  }

}


#' Matching Gene IDs
#'
#' Makes sure the input_data has the same gene IDs in the same order as
#' the ADAGE model. It first converts gene IDs from the input_data to PAO1 locus
#' tags. Then it re-orders input's rows according to gene order in ADAGE and fills
#' in zero values for genes used in the ADAGE model but missed in the
#' input_data.
#'
#' @param input_data a data.frame containing geneIDs in the first column and each
#' sample's gene expression values from the second column.
#' @param ref_IDs a vector storing the reference geneIDs in the right order
#' @return A data.frame containing the input_data's expression values after
#' converting gene IDs, sorting gene orders, and filling in missing genes.
match_IDs <- function(input_data, ref_IDs){

  if (!check_input(input_data)) {
    stop("The input data should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric
         expression values for each sample.")
  }

  # convert the gene IDs used in the input_data to PAO1 locus tags
  converted_geneIDs <- sapply(as.data.frame(input_data)[, 1],
                              function(x) to_LocusTag(input_ID = x,
                                                      ref_IDs = ref_IDs))

  # create a match index between input IDs and reference IDs
  match_index <- match(ref_IDs, converted_geneIDs)

  # NAs are created if the input does not contain some reference IDs
  na_index <- which(is.na(match_index))

  # print a warning if some IDs are missed in the input
  if (length(na_index) > 0) {
    na_geneID <- ref_IDs[na_index]
    warning(paste("ADAGE gene features", paste(na_geneID, collapse = ","),
                  "are not found in the input!"))
  }

  # re-order the input_data
  IDmapped <- input_data[match_index, -1]

  # assign reference ID to the input_data
  IDmapped <- data.frame(geneID = ref_IDs, IDmapped, stringsAsFactors = FALSE,
                         check.names = FALSE)

  return(IDmapped)
}


#' Missing values imputation
#'
#' Imputes and fills missing values using k-nearest neighbor method. The 5
#' nearest neighbors of missing genes are calculated mainly using the reference
#' data and then are used to fill in missing values in the input_data. It uses
#' the impute.knn function from the impute bioconductor package.
#'
#' @param input_data a data.frame with gene IDs in the first column and
#' expression values from the second column.
#' @param ref_data a data.frame storing gene expression values to be used
#' as a reference dataset for calculating nearest neighbors. Should have more
#' samples than the input_data. The first columns of input_data and ref_data that
#' specify gene IDs should exactly be the same.
#' @return the input_data data.frame with missing values being filled
#' @seealso \code{\link[impute]{impute.knn}}
impute_miss_values <- function(input_data, ref_data){

  if (!check_input(input_data)) {
    stop("The input data should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric
         expression values for each sample.")
  }

  if (!check_input(ref_data)) {
    stop("The reference data should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric
         expression values for each sample.")
  }

  # make sure each row in the input_data and reference data represents the
  # same gene (the first columns are the same)
  if (!(nrow(input_data) == nrow(ref_data)) |
      !all(input_data[, 1] == ref_data[, 1])) {
    stop("The gene identifiers from the input_data and reference data should be
         the same!")
  }

  # combine the input_data and the reference data
  combined_data <- as.matrix(cbind(input_data[, -1], ref_data[, -1]))
  # perform knn imputation using the reference data to find nearest neighbors
  impute_result <- impute::impute.knn(combined_data, k = 5)
  # extract the imputed input_data
  imputed_data <- impute_result$data[, 1:ncol(input_data[, -1])]
  # add geneID column in front
  imputed_data <- data.frame(geneID = input_data[, 1], imputed_data,
                             stringsAsFactors = FALSE, check.names = FALSE)
  return(imputed_data)
}


#' Quantile normalization
#'
#' Performs quantile normalization on the input_data. If using reference,
#' quantile normalization is done using the quantile distribution derived
#' from the reference data.
#'
#' @param input_data a data.frame with gene IDs in the first column and
#' expression values from the second column.
#' @param use_ref a logical value indicating whether the normalization should
#' be done based on the reference data (default: FALSE).
#' @param ref_data a data.frame storing gene expression values to be used
#' as a reference dataset for quantile normalization. The first columns of
#' input_data and ref_data that specify gene IDs should exactly be the same.
#' @return the input_data.frame after quantile normalization.
#' @seealso \code{\link[preprocessCore]{normalize.quantiles}},
#' \code{\link[preprocessCore]{normalize.quantiles.use.target}},
#' \code{\link[preprocessCore]{normalize.quantiles.determine.target}}
quantile_norm <- function(input_data, use_ref = FALSE, ref_data){

  if (!check_input(input_data)) {
    stop("The input data should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric
         expression values for each sample.")
  }

  if (!check_input(ref_data)) {
    stop("The reference data should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric
         expression values for each sample.")
  }

  # make sure each row in the input_data and reference data represents the
  # same gene (the first columns are the same)
  if (!(nrow(input_data) == nrow(ref_data)) |
      !all(input_data[, 1] == ref_data[, 1])) {
    stop("The gene identifiers from the input_data and reference data should be
         the same!")
  }

  if (use_ref) {
    # determine quantile distribution using the reference data
    ref_dist <- preprocessCore::normalize.quantiles.determine.target(
      as.matrix(ref_data[, -1]))
    # normalize the input_data with the derived quantile distribution
    normed_data <- preprocessCore::normalize.quantiles.use.target(
      as.matrix(input_data[, -1]), ref_dist)
  } else {
    # directly normalize the input_data
    normed_data <- preprocessCore::normalize.quantiles(
      as.matrix(input_data[, -1]))
  }

  # add geneID column in front
  normed_data <- data.frame(geneID = input_data[, 1], normed_data,
                            stringsAsFactors = FALSE, check.names = FALSE)
  colnames(normed_data) <- colnames(input_data)

  return(normed_data)
}


#' RNAseq data normalization with TDM
#'
#' Performs Training Distribution Matching (TDM) on the input RNAseq data
#' to normalize RNAseq expression values to a comparable range of the reference
#' microarray data.
#'
#' @param input_data a data.frame with gene IDs in the first column and
#' expression values from the second column.
#' @param ref_data a data.frame storing gene expression values of a
#' compendium built from microarray data. It is used as the reference data in
#' TDM.
#' @return a data.frame storing TDM normalized gene expression values from the
#' input_data.
#' @seealso \url{https://github.com/greenelab/TDM}
TDM_RNAseq <- function(input_data, ref_data){

  if (!check_input(input_data)) {
    stop("The input data should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric
         expression values for each sample.")
  }

  if (!check_input(ref_data)) {
    stop("The reference data should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric
         expression values for each sample.")
  }

  # make sure each row in the input_data and reference data represents the
  # same gene (the first columns are the same)
  if (!(nrow(input_data) == nrow(ref_data)) |
      !all(input_data[, 1] == ref_data[, 1])) {
    stop("The gene identifiers from the input_data and reference data should be
         the same!")
  }

  # TDM require the first column to be named as "gene" and use data.table
  # data structure
  colnames(input_data)[1] <- "gene"
  colnames(ref_data)[1] <- "gene"
  input_data <- data.table::data.table(input_data)
  data.table::setkey(input_data, gene)
  ref_data <- data.table::data.table(ref_data)
  data.table::setkey(ref_data, gene)

  # perform TDM
  data_tdm <- TDM::tdm_transform(input_data, ref_data)

  # convert data.table back to data.frame
  data.table::setDF(data_tdm)
  colnames(data_tdm)[1] <- "geneID"
  data_tdm$geneID <- as.character(data_tdm$geneID)
  data_tdm[, -1] <- as.numeric(as.matrix(data_tdm[, -1]))

  return(data_tdm)

}


#' Zero-one normalization
#'
#' Normalizes gene expression values in the input_data to be between 0 and 1.
#' Normalization is done per gene (row-wise) through
#' subtracting row minimum from each value and then being divided by row range.
#' If using reference, the normalization is done using a gene's expression
#' minimum and range in the reference data.
#'
#' @param input_data a data.frame with gene IDs in the first column and
#' expression values from the second column.
#' @param use_ref a logical value indicating whether the normalization should
#' be done based on the reference data (default: FALSE).
#' @param ref_data a data.frame storing gene expression values to be used
#' as a reference dataset for zero-one normalization. The first columns in
#' input_data and ref_data that specify gene IDs should be exactly the same.
#' @return a data.frame storing zero-one normalized gene expression values from
#' the input_data.
#' @export
zeroone_norm <- function(input_data, use_ref = FALSE, ref_data) {

  if (!check_input(input_data)) {
    stop("The input data should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric
         expression values for each sample.")
  }

  if (!check_input(ref_data)) {
    stop("The reference data should be a data.frame with the first column as a
         character of gene IDs and the rest of the columns storing numeric
         expression values for each sample.")
  }

  # make sure each row in the input_data and reference data represents the
  # same gene (the first columns are the same)
  if (!(nrow(input_data) == nrow(ref_data)) |
      !all(input_data[, 1] == ref_data[, 1])) {
    stop("The gene identifiers from the input_data and reference data should be
         the same!")
  }

  if (use_ref) {

    # get gene expression ranges and minima from the reference data
    ref_range <- apply(ref_data[, -1], 1, function(x) diff(range(x)))
    ref_min <- apply(ref_data[, -1], 1, min)

    # perform zero-one normalization using the ranges and minima from the
    # reference data
    zeroone_normed <- (input_data[, -1] - ref_min) / ref_range

    # bound the value to be between 0 and 1
    zeroone_normed[zeroone_normed > 1] <- 1
    zeroone_normed[zeroone_normed < 0] <- 0

  } else {

    # perform zero-one normalization per gene directly on the input_data
    zeroone_normed <- t(apply(input_data[, -1], 1,
                              function(x) (x - min(x)) / diff(range(x))))

  }

  # build the output data.frame
  zeroone_normed <- data.frame(geneID = as.data.frame(ref_data)[, 1],
                               zeroone_normed, stringsAsFactors = FALSE,
                               check.names = FALSE)

  return(zeroone_normed)
}


#' Checking input format
#'
#' Checks whether the input is a data.frame (or a tibble) with its first
#' column being character and the rest columns being numeric.
#'
#' @param input_data the input to check
#' @return TRUE if the input_data.frame meets the requirements, otherwise FALSE.
check_input <- function(input_data){

  # check whether the first column is character and the rest columns
  # are numeric.

  if (is.data.frame(input_data)) {

    # use as.data.frame if input_data is a tibble
    if (is.character(as.data.frame(input_data)[, 1]) &
        all(sapply(input_data[, -1], is.numeric))) {

      return(TRUE)

    }
  }

  return(FALSE)

}
