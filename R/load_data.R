#' Data loading and preprocessing
#'
#' Currently supports three types of inputs:
#' a zip file containing microarray CEL files in a dataset;
#' a folder containing microarray CEL files in a dataset;
#' a tab-delimited txt file storing processed gene expression values in a
#' dataset (gene identifiers in the first column and then one sample per column).
#'
#' @param input file path to the input file or folder.
#' @param isProcessed a logical value indicating whether the input data has
#' already been processed (default: FALSE).
#' @param isRNAseq a logical value indicating whether the input data is RNAseq
#' data. It determines whether TDM is applied to normalize the expression
#' values (default: FALSE).
#' @return A data frame containing the processed gene expression values ready
#' for ADAGE analysis.
#' @examples
#' load_data(filepath, isProcessed = FALSE, isRNAseq = FALSE)
#' @importFrom readr read_tsv
load_dataset <- function(input, isProcessed = FALSE, isRNAseq = FALSE){

  # quantile normalize probes if the input is raw CEL formot, directly load data
  # if the input is processed microarray or RNAseq data
  if (!isProcessed) {
    # if the input is a zip file, unzip into the same folder and use cel files
    # in the folder. If the input is a folder, directly use cel files in the
    # folder

    if (endsWith(input, ".zip")) {
      # unzip files
      unzip(input, exdir = dirname(input))
      cel_folder <- file.path(dirname(input), gsub(".zip", "", basename(input)))

    } else if (R.utils::isDirectory(input)) {
      cel_folder <- input
    }

    # normalize each probe
    data <- process_celfiles(cel_folder = cel_folder, use_ref = TRUE,
                             quantile_ref = probedistribution)

  } else {

    # read in the processed data
    data <- read_tsv(input)
    colnames(data)[1] <- "geneID"

  }

  # transform gene features in the input data to gene features used in ADAGE model
  data <- match_IDs(input_data = data, ref_IDs = eADAGEmodel$geneID)

  # perform TDM transformation if the input is RNAseq data
  if (isRNAseq) {
    data <- TDM_RNAseq(input_data = data, ref_data = compendium)
  }

  # zero-one normalization
  data <- zeroone_norm(input_data = data, use_ref = TRUE, ref_data = compendium)

  return(data)
}


#' CEL files processing
#'
#' Processes microarray data in cel format. Only cel files measured on the
#' "Pae_G1a" platform will be processed. If using reference, the quantile
#' normalization of probes is based on a reference quantile
#' distribution.
#'
#' @param cel_folder file path to the folder storing cel files.
#' @param use_ref a logical value indicating whether the probe normalization
#' should be done using a reference quantile distribution (default: TRUE).
#' @param quantile_ref a vector storing the reference quantile distribution
#' (default: the quantile distribution of probes used in
#' normalzing the P.a. gene expression compendium) Since the input microarray
#' data need to be normalzied to the processed compendium, the default
#' should not be changed in most cases.
#' @return a data frame containing normalized gene expression values with geneID
#' in the first column and then each cel file in one column.
process_celfiles <- function(cel_folder, use_ref = TRUE,
                             quantile_ref = probedistribution){

  # read in all cel files in the folder
  celfiles <- affy::list.celfiles(cel_folder)

  # get the platform types
  ptype <- sapply(celfiles, function(f) {
    affyio::read.celfile.header(paste(cel_folder, f, sep = "/"))[1]
  })

  # only process cel files on "Pae_G1a" platform
  pfiles <- paste(cel_folder, subset(celfiles, ptype == "Pae_G1a"), sep = "/")

  # read in cel files and produce an AffyBatch object
  affybatch <- affy::ReadAffy(filenames = pfiles)


  if (use_ref){

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

    # directly use rma function to normalize within the input data
    expression <- affy::rma(affybatch)

  }

  # build the final gene expression data frame
  colnames(expression) <- celfiles
  expression <- data.frame(geneID = rownames(expression), expression,
                           stringsAsFactors = FALSE)

  return(expression)

}


#' Gene ID convertion
#'
#' Converts the input gene ID to PAO1 locus tag that are used in ADAGE model.
#' It currently can convert gene IDs on AFFY chip, PAO1 gene symbols, and
#' PAO1 gene orthologs in other P.a. strains. It returns NA if the input is not
#' recognized as one of above.
#'
#' @param input_ID the input gene ID (character).
#' @return the corresponding PAO1 locus tag ("PAXXXX") for the input gene or NA
#' if the input gene is not recognized.
to_LocusTag <- function(input_ID) {

  if (input_ID %in% eADAGEmodel$geneID) {

    # do nothing if the input ID is already used in eADAGE model
    return(input_ID)

  } else if (startsWith(input_ID, "ig") | startsWith(input_ID, "Pae") |
      startsWith(input_ID, "AFFY")) {

    # gene IDs start with "ig", "Pae", or "AFFY" are controls on the AFFY chip
    return("control")

  } else if (endsWith(input_ID, "_at")) {

    # gene IDs that end with "_at" but are not controls are AFFY IDs in the
    # format "PAXXXX_symbol", only the first part before "_" will be perserved
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


#'Gene IDs matching
#'
#'Makes sure the input data having the same gene IDs in the same order as
#'the ADAGE model. It first converts gene IDs from the input data to PAO1 locus
#'tags. Then it re-orders input's rows according to gene order in ADAGE and fills
#'in zero values for genes used in the ADAGE model but missed in the
#'input data.
#'
#'@param input_data a data frame containing geneIDs in the first column and each
#'sample's gene expression values from the second column.
#'@param ref_IDs a vector storing the reference geneIDs in the right order
#'(default: ordered PAO1 locus tags used in the ADAGE model).
#'@return A data frame containing the input_data's expression values after
#'converting gene IDs, sorting gene orders, and filling in missing genes.
match_IDs <- function(input_data, ref_IDs = eADAGEmodel$geneID){

  # convert the gene IDs used in the input data to PAO1 locus tags
  converted_geneIDs <- sapply(input_data$geneID, function(x) to_LocusTag(x))

  # create a match index between input IDs and reference IDs
  match_index <- match(ref_IDs, converted_geneIDs)

  # NAs are created if the input does not contain some reference IDs
  na_index <- which(is.na(match_index))

  # print a warning if some IDs are missed in the input
  if (length(na_index) > 0){
    na_geneID <- ref_IDs[na_index]
    warning(paste("ADAGE gene features", paste(na_geneID, collapse = ","),
                  "are not found in the input!",
                  "Their expression values will be set to 0."))
  }

  # re-order the input data
  IDmapped <- input_data[match_index, ]
  # set NA to 0
  IDmapped[is.na(IDmapped)] <- 0
  # assign reference ID to the input data
  IDmapped$geneID <- ref_IDs

  return(IDmapped)
}


#'RNAseq data normalization with TDM
#'
#'Performs Training Distribution Matching (TDM) on the input RNAseq data
#'to normalize RNAseq data to a comparable range of the reference microarray
#'compendium.
#'
#'@param input_data A data frame with gene IDs in the first column and
#'expression values from the second column.
#'@param ref_data A data frame storing gene expression values of a
#'compendium built from microarray data. It is used as the reference data in
#'TDM.
#'@return A data frame storing TDM normalized gene expression values from the
#'input data.
#'@seealso \url{https://github.com/greenelab/TDM}
TDM_RNAseq <- function(input_data, ref_data = compendium){

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
#'Normalizes the input gene expression values to be between 0 and 1.
#'Normalization is done per gene (row) through
#'substracting its minimun from each value and then being divided by its range.
#'If using reference, the normalization is done using the gene's expression
#'mininum and range in the reference data.
#'
#'@param input_data A data frame with gene IDs in the first column and
#'expression values from the second column.
#'@param use_ref A logical value indicating whether the normalization should
#'be done based on the reference data.
#'@param ref_data A data frame storing gene expression values to be used
#'as a reference dataset for zero-one normalization (default: the P.a. gene
#'expression compendium). The first columns in input_data and ref_data that
#'specify gene IDs should be exactly the same.
#'@return A data frame storing zero-one normalized gene expression values from
#'the input data.
zeroone_norm <- function(input_data, use_ref = FALSE, ref_data = compendium) {

  # make sure each row in the input data and reference data represents the
  # same gene (the first columns are the same)
  if (!(nrow(input_data) == nrow(ref_data)) |
      !all(input_data[, 1] == ref_data[, 1])) {
    stop("The gene identifiers from the input data and reference data should be
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

    # perform zero-one normalization per gene directly on the input data
    zeroone_normed <- apply(input_data[, -1], 1,
                            function(x) (x - min(x)) / diff(range(x)))

  }

  # build the output data frame
  zeroone_normed <- data.frame(geneID = input_data$geneID, zeroone_normed,
                               stringsAsFactors = FALSE)

  return(zeroone_normed)
}

