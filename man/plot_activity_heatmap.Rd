% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cal_activity.R
\name{plot_activity_heatmap}
\alias{plot_activity_heatmap}
\title{Signature activity heatmap}
\usage{
plot_activity_heatmap(activity, signatures = NULL, fix_color_range = TRUE,
  is_pathway = FALSE)
}
\arguments{
\item{activity}{a data.frame that stores the signature activities for each
sample in a dataset. The first column is signature name and activity values
start from the second column.}

\item{signatures}{a character vector specifying which signatures to include
in the heatmap (default: NULL, all signatures will be included).}

\item{fix_color_range}{logical. If TRUE, fix the heatmap color
range to the maximum activity range of the input activity. If FALSE, heatmap
color range is determined by the activity ranges of the selected signatures.
(default: TRUE)}

\item{is_pathway}{logical, indicates whether the input activity stores
pathway activity. If TRUE, figure margin is set to a larger value to better
display pathway names, which are much longer than signature names.}
}
\description{
Plots a heatmap showing signature activities in a dataset.
}

