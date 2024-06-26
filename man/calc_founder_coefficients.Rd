% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc-founder-coefficients.R
\name{calc_founder_coefficients}
\alias{calc_founder_coefficients}
\title{Get the founder coefficients.}
\usage{
calc_founder_coefficients(
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

\item{blup}{`TRUE` to perform BLUP.}

\item{center}{`TRUE` to center the data.}

\item{cores}{The number of cores to use (0 = ALL).}
}
\value{
A named `list` with each element being a tibble with the following
}
\description{
Get the founder coefficients.
}
