#' eADAGE model
#'
#' A 300-node eADAGE model built for biological pathway analysis for
#' P. aeruginosa
#'
#' @format A data frame with 5549 rows and 301 variables:
#' \describe{
#'   \item{geneID}{P.a. PAO1 gene identifiers}
#'   \item{Node1}{Each gene's weight to Node1}
#'   ...
#' }
#' @source \url{http://biorxiv.org/content/early/2016/10/03/078659}
"eADAGEmodel"

#' probe distribution
#'
#' The quantile distribution of microarray probes used in building the
#' P. aeruginosa compendium. It will be used as a reference distribution to
#' quantile normalize new datasets.
#'
#' @format A numeric vector:
"probedistribution"

#' gene expression range
#'
#' The expression range of each gene in the P. aeruginosa compendium.
#'
#' @format A data frame with 5549 rows and 3 variables:
#' \describe{
#'  \item{geneID}{P.a. PAO1 gene identifiers}
#'  \item{max_express}{Each gene's maximum expression value across the compendium}
#'  \item{min_express}{Each gene's minimum expression value across the compendium}
#' }
"expressionrange"

#' gene information
#'
#' Pseudomonas aeruginosa gene information curated by NCBI. It is used for
#' mapping PA numbers (LocusTag) to gene Symbol.
#'
#' @format A data frame with 5698 rows and 15 variables:
#' \describe{
#'  \item{#taxid}{Taxonomic id of P. aeruginosa}
#'  \item{GeneID}{Unique gene identifiers}
#'  \item{Symbol}{Gene symbols or gene names}
#'  \item{LocusTag}{Locus tags or PA numbers}
#'  ...
#' }
#' @source \url{ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Archaea_Bacteria/Pseudomonas_aeruginosa_PAO1.gene_info.gz}
"geneinfo"
