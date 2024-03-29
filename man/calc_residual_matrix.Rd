% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-correlation.R
\name{calc_residual_matrix}
\alias{calc_residual_matrix}
\title{Calculate the residual matrix.}
\usage{
calc_residual_matrix(
  variable_matrix,
  adjust_matrix,
  variables_interest,
  variables_compare,
  use_qr = TRUE
)
}
\arguments{
\item{variable_matrix}{The data \code{matrix} for first set.}

\item{adjust_matrix}{The data \code{matrix} for the second set.}

\item{variables_interest}{\code{list} of variables of interest.}

\item{variables_compare}{\code{list} of variables to compare.}

\item{use_qr}{\code{TRUE} to use qr decomposition.}
}
\value{
The residual \code{matrix}.
}
\description{
Calculate the residual matrix.
}
