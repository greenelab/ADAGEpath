% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_signature.R
\name{build_gene_geneset_map}
\alias{build_gene_geneset_map}
\title{Gene-geneset map}
\usage{
build_gene_geneset_map(genesets)
}
\arguments{
\item{genesets}{a named list with each element storing genes in one
geneset.}
}
\value{
a data.frame with the first column specifying geneID and the second
column specifying the names of genesets that a gene is in.
}
\description{
Map genes to the genesets they are in.
}

