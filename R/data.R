#' eADAGE model
#'
#' The 300-node eADAGE model used in the paper "Unsupervised extraction of
#' functional gene expression signatures in the bacterial pathogen Pseudomonas
#' aeruginosa with eADAGE"
#'
#' @format A data frame (tibble) with 5549 rows and 301 variables. The first
#' column is geneID that specifies gene identifiers. Starting from
#' the second column, each column stores numeric weight values of one node
#' @source \url{https://doi.org/10.1101/078659}
"eADAGEmodel"


#' Pseudomonas aeruginosa gene expression compendium
#'
#' The compendium contains all P.a. microarray datasets measured on the
#' [Pae_G1a] Affymetric Pseudomonas aeruginosa array
#' platform from the ArrayExpress database as of 31 July 2015. The expression
#' values in the compendium have been background corrected and quantile
#' normalized.
#'
#' @format A data frame (tibble) with 5549 rows and 1052 variables. The first
#' column is geneID that specifies gene identifiers. Starting from
#' the second column, each column stores numeric expression values of one sample.
#' @source \url{https://doi.org/10.1101/078659}
"PAcompendium"


#' Probe distribution
#'
#' The quantile distribution of microarray probes used in building the preloaded
#' PAcompendium data object. It will be used as a reference distribution to
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
#' @format A data.frame (tibble) with 5698 rows and 15 variables:
#' \describe{
#'   \item{\samp{#tax_id}}{a numeric vector, taxonomic id of P. aeruginosa}
#'   \item{\code{GeneID}}{a numeric vector, unique gene identifiers}
#'   \item{\code{Symbol}}{a character vector, gene symbols or gene names}
#'   \item{\code{LocusTag}}{a character vector, locus tags or PA numbers}
#'   \item{\code{Synonyms}}{a numeric vector}
#'   \item{\code{dbXrefs}}{a character vector}
#'  \item{\code{chromosome}}{a numeric vector}
#'   \item{\code{map_location}}{a numeric vector}
#'   \item{\code{description}}{a character vector}
#'   \item{\code{type_of_gene}}{a character vector}
#'   \item{\code{Symbol_from_nomenclature_authority}}{a numeric vector}
#'   \item{\code{Full_name_from_nomenclature_authority}}{a numeric vector}
#'   \item{\code{Nomenclature_status}}{a numeric vector}
#'   \item{\code{Other_designations}}{a numeric vector}
#'   \item{\code{Modification_date}}{a numeric vector}
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
#' @format A data.frame (tibble) with 392139 rows and 8 variables:
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


#' Experiments and samples in Pseudomonas aeruginosa gene expression compendium
#'
#' ArrayExpress experiment assession numbers and their associated sample IDs
#' included in the preloaded PAcompendium data object.
#' There are duplicates in sample IDs because one sample can be included in
#' multiple experiments.
#' @format A data.frame with 1118 rows and 2 variables:
#' \describe{
#'  \item{Experiment}{ArrayExpress experiment assession numbers, only contain
#'  experiments in the preloaded PAcompendium data object.}
#'  \item{Sample}{Sample IDs in "xxx.CEL" format}
#' }
"experimentID"


#' Pseudomonas aeruginosa PAO1 operons
#'
#' Pseudomonas aeruginosa PAO1 operon predictions obtained from the DOOR (
#' Database of prOkaryotic OpeRons) database on Dec. 12 2016.
#'
#' @format A list with each element being a character vector that stores genes
#' in one operon. Locus tag (PAXXXX) is used as gene ID.
#' @source \url{http://csbl.bmb.uga.edu/DOOR/index.php}
"operons"
