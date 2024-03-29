% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-founder-coefficients.R
\name{get_founder_coefficients}
\alias{get_founder_coefficients}
\title{Get the founder coefficients.}
\usage{
get_founder_coefficients(
  dataset,
  id,
  chrom,
  intcovar = NULL,
  blup = FALSE,
  center = TRUE,
  cores = 0
)
}
\arguments{
\item{dataset}{The dataset object.}

\item{id}{The unique id in the dataset.}

\item{chrom}{The chromosome.}

\item{intcovar}{The interactive covariate.}

\item{blup}{\code{TRUE} to perform BLUP.}

\item{center}{\code{TRUE} to center the data.}

\item{cores}{The number of cores to use (0 = ALL).}
}
\value{
A named \code{list} with each element being a tibble with the following
}
\description{
Get the founder coefficients.
}
