% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-lod-scan.R
\name{get_lod_scan}
\alias{get_lod_scan}
\title{Perform a LOD scan.}
\usage{
get_lod_scan(
  dataset,
  id,
  intcovar = NULL,
  cores = 0,
  filter_threshold = 6,
  filter_peak_drop = Inf,
  filter_thresholdX = NULL,
  filter_peak_dropX = NULL,
  scan1_output = FALSE
)
}
\arguments{
\item{dataset}{The dataset object.}

\item{id}{The unique id in the dataset.}

\item{intcovar}{The interactive covariate.}

\item{cores}{Number of cores to use (0 = ALL).}

\item{filter_threshold}{If set, qtl2::find_peaks is used.}

\item{filter_peak_drop}{If set, qtl2::find_peaks is used.}

\item{filter_thresholdX}{If set, qtl2::find_peaks is used.}

\item{filter_peak_dropX}{If set, qtl2::find_peaks is used.}

\item{scan1_output}{If \code{TRUE}, original \code{qtl2::scan1} data is included.}
}
\value{
a \code{list} with the following elements:
\itemize{
\item lod_peaks - \code{tibble} of LOD peaks
\item lod_scores - \code{tibble} with the following columns: id, chr, pos, lod
\item scan1 - \code{qtl2::scan1} output
}
}
\description{
Perform a LOD scan.
}
