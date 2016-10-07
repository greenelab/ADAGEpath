load_compendium <- function(compendium_file, normed = TRUE) {

}

load_model <- function(model_file, geneID_file){

}

convert_geneID <- function(IDmap_file, key){

}

load_dataset <- function(dataset_folder, selected_samples, phenotypes,
                         quantile_ref, output_file){
  qnormed <- normalize_celfiles(dataset_folder, quantile_ref)
  qnormed_PAonly <- remove_nonPAIDs(qnormed)
  qnormed_PAonly_01normed <- zeroone_norm(qnormed_PAonly)

}


#' @importFrom readr read_delim
normalize_celfiles <- function(cel_folder, quantile_ref){
  # read in the reference quantile distribution
  quantile_ref <- read_delim(quantile_ref, col_names = FALSE, delim = " ")
  quantile_ref <- quantile_ref[[1]]

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
  PMmat_sumed <- cbind(gene = rownames(PMmat_sumed), PMmat_sumed)

  return(PMmat_sumed)

}


remove_nonPAIDs <- function() {

}

zeroone_norm <- function() {

}
