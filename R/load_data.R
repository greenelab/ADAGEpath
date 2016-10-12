
load_dataset <- function(input){
  if (endsWith(input, ".zip")) {
    # unzip files
    unzip(input, exdir = dirname(input))
    cel_folder <- file.path(dirname(input), gsub(".zip", "", basename(cel_zip)))
  } else if (R.utils::isDirectory(input)) {
    cel_folder <- input
  }
  qnormed <- normalize_celfiles(cel_folder = cel_folder,
                                quantile_ref = probedistribution)
  geneIDs <- sapply(qnormed$geneID, function(x) to_LocusTag(x))
  qnormed_IDmapped <- qnormed[match(eADAGEmodel$geneID, geneIDs), ]
  qnormed_IDmapped$geneID <- eADAGEmodel$geneID
  return(qnormed_IDmapped)
}


#' @importFrom readr read_delim
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

TDM_RNAseq <- function(){

}

zeroone_norm <- function() {

}


