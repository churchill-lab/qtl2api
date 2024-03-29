% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-lod-peaks.R
\name{get_lod_peaks_for_annot}
\alias{get_lod_peaks_for_annot}
\title{Get all the peaks for an id.}
\usage{
get_lod_peaks_for_annot(
  dataset,
  id,
  threshold = 6,
  peakdrop = 2,
  thresholdX = 6,
  peakdropX = 2,
  n_cores = 0
)
}
\arguments{
\item{dataset}{the dataset object}

\item{id}{the unique id in the dataset}

\item{threshold}{if set, qtl2::find_peaks is used}

\item{peakdrop}{if set, qtl2::find_peaks is used}

\item{thresholdX}{if set, qtl2::find_peaks is used}

\item{peakdropX}{if set, qtl2::find_peaks is used}

\item{n_cores}{number of cores to use (0=ALL)}
}
\value{
a tibble of the peaks.
}
\description{
Get all the peaks for an id.
}
