% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_data.R
\name{match_IDs}
\alias{match_IDs}
\title{Matching Gene IDs}
\usage{
match_IDs(input_data, ref_IDs)
}
\arguments{
\item{input_data}{a data.frame containing geneIDs in the first column and each
sample's gene expression values from the second column.}

\item{ref_IDs}{a vector storing the reference geneIDs in the right order}
}
\value{
A data.frame containing the input_data's expression values after
converting gene IDs, sorting gene orders, and filling in missing genes.
}
\description{
Makes sure the input_data has the same gene IDs in the same order as
the ADAGE model. It first converts gene IDs from the input_data to PAO1 locus
tags. Then it re-orders input's rows according to gene order in ADAGE and fills
in zero values for genes used in the ADAGE model but missed in the
input_data.
}

