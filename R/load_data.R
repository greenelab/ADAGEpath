
#' @importFrom readr read_tsv
load_dataset <- function(input, isProcessed = FALSE, isRNAseq = FALSE){

  # quantile normalize probes if the input is raw CEL files, directly load data
  # if the input is processed microarray or RNAseq data
  if (!isProcessed) {
    # if the input is a zip file
    if (endsWith(input, ".zip")) {
      # unzip files
      unzip(input, exdir = dirname(input))
      cel_folder <- file.path(dirname(input), gsub(".zip", "", basename(input)))
    } else if (R.utils::isDirectory(input)) {
      cel_folder <- input
    }
    data <- normalize_celfiles(cel_folder = cel_folder,
                               quantile_ref = probedistribution)
  } else {
    # read in the processed data
    data <- read_tsv(input)
    colnames(data)[1] <- "geneID"
  }

  # transform gene IDs in the input data to gene IDs used in ADAGE model
  data <- match_IDs(input_data = data, ref_IDs = eADAGEmodel$geneID)

  # perform TDM transformation if the input is RNAseq data
  if (isRNAseq) {
    data <- TDM_RNAseq(input_data = data, ref_data = compendium)
  }

  # zero-one normalization
  data <- zeroone_norm(input_data = data, use_ref = TRUE, ref_data = compendium)
  return(data)
}


normalize_celfiles <- function(cel_folder, quantile_ref = probedistribution){

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
  PMmat_sumed <- preprocessCore::subColSummarizeMedianpolishLog(
    PMmat_normed, group.labels = probe_list)

  # build the final gene expression matrix
  colnames(PMmat_sumed) <- celfiles
  PMmat_sumed <- data.frame(geneID = rownames(PMmat_sumed), PMmat_sumed,
                            stringsAsFactors = FALSE)

  return(PMmat_sumed)

}

to_LocusTag <- function(input_ID) {
  if (startsWith(input_ID, "ig") | startsWith(input_ID, "Pae") |
      startsWith(input_ID, "AFFY")) {
    output_ID <- "control"
  } else if (endsWith(input_ID, "_at")) {
    output_ID <- unlist(strsplit(input_ID, "_"))[1]
  } else if (input_ID %in% geneinfo$Symbol) {
    output_ID <- geneinfo$LocusTag[geneinfo$Symbol == input_ID]
  } else if (input_ID %in% geneinfo$LocusTag) {
    output_ID <- input_ID
  } else if (input_ID %in% PAO1orthologs$`Locus Tag (Hit)`) {
    output_ID <- PAO1orthologs$`Locus Tag (Query)`[
      PAO1orthologs$`Locus Tag (Hit)` == input_ID]
  } else {
    warning("Gene ID not found!")
    output_ID <- NA
  }
  return(output_ID)
}

match_IDs <- function(input_data, ref_IDs = eADAGEmodel$geneID){
  converted_geneIDs <- sapply(input_data$geneID, function(x) to_LocusTag(x))
  match_index <- match(ref_IDs, converted_geneIDs)

  na_index <- which(is.na(match_index))
  if (length(na_index) > 0){
    na_geneID <- ref_IDs[na_index]
    warning(paste(paste(na_geneID, collapse = ","), "are not found in the input!",
                  "Their expression values will be set to 0."))
  }

  IDmapped <- input_data[match_index, ]
  # set NA to 0
  IDmapped[is.na(IDmapped)] <- 0
  IDmapped$geneID <- ref_IDs
  return(IDmapped)
}


TDM_RNAseq <- function(input_data, ref_data = compendium){

  # install and load the TDM package from github
  devtools::install_github("greenelab/TDM")
  library(TDM)

  input_data <- data.table::data.table(input_data)
  data.table::setkey(input_data, geneID)
  ref_data <- data.table::data.table(ref_data)
  data.table::setkey(ref_data, geneID)

  data_tdm <- tdm_transform(input_data, ref_data)
  return(data_tdm)

}

zeroone_norm <- function(input_data, use_ref = FALSE, ref_data = compendium) {

  # make sure each row in the input data and reference data represents the
  # same gene
  if (input_data[, 1] != ref_data[, 1]) {
    stop("The gene identifiers from the input data and reference data should be
         the same!")
  }

  if (use_ref) {
    ref_range <- apply(ref_data[, -1], 1, function(x) diff(range(x)))
    ref_min <- apply(ref_data[, -1], 1, min)
    zeroone_normed <- (input_data[, -1] - ref_min) / ref_range
    # bound the value to be between 0 and 1
    zeroone_normed[zeroone_normed > 1] <- 1
    zeroone_normed[zeroone_normed < 0] <- 0
  } else {
    zeroone_normed <- apply(input_data[, -1], 1,
                            function(x) (x - min(x)) / diff(range(x)))
  }

  zeroone_normed <- data.frame(geneID = input_data$geneID, zeroone_normed,
                               stringsAsFactors = FALSE)
  return(zeroone_normed)
}


