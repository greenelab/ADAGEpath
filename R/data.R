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


#' Pseudomonas aeruginosa gene expression compendium
#'
#' The compendium contains all P.a. microarray datasets measured on GPL84
#' platform from the ArrayExpress database as of 31 July 2015. The expression
#' values in the compendium have been background corrected and quantile
#' normalized.
#'
#' @format A data frame with 5549 rows and 1052 variables. The first column
#' specifies gene identifiers and samples start from the second column.:
#' \describe{
#'   \item{geneID}{P.a. PAO1 gene identifiers}
#'   ...
#' }
#' @source \url{http://biorxiv.org/content/early/2016/10/03/078659}
"compendium"


#' Probe distribution
#'
#' The quantile distribution of microarray probes used in building the
#' P. aeruginosa compendium. It will be used as a reference distribution to
#' quantile normalize new datasets.
#'
#' @format A numeric vector:
"probedistribution"


#' Gene information
#'
#' Pseudomonas aeruginosa gene information curated by NCBI.
#' It was downloaded from NCBI ftp on Oct. 7 2016. It is used for
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


#' Pseudomonas aeruginosa PAO1 strain orthologs
#'
#' Pseudomonas aeruginosa PAO1 strain ortholog predictions obtained from
#' http://pseudomonas.com on Oct. 10 2016.
#' It uses genes in PAO1 strain as query and finds
#' orthologs hit in other P.a. strains.
#'
#' @format A data frame with 392139 rows and 8 variables:
#' \describe{
#'  \item{Strain (Query)}{Pseudomonas aeruginosa PAO1 (Reference)}
#'  \item{Locus Tag (Query)}{Locus tag for PAO1}
#'  \item{Description (Query)}{Description of the PAO1 gene}
#'  \item{Strain (Hit)}{Ortholog hit strain}
#'  \item{Locus Tag (Hit)}{Locus tag of the hit ortholog}
#'  \item{Description (Hit)}{Description of the hit ortholog}
#'  \item{Percent Identity}{Percentage of sequence identity}
#'  \item{Alignment Length}{Alignment length}
#'  ...
#' }
#' @source \url{http://pseudomonas.com/downloads/pseudomonas/pgd_r_16_1/Pseudomonas_aeruginosa_PAO1_107/Pseudomonas_aeruginosa_PAO1_107_orthologs.txt}
"PAO1orthologs"

