# #############################################################################
#
# Utility functions
#
# #############################################################################

`%not in%` <- function(x, table) match(x, table, nomatch = 0L) == 0L


#' Get the version on qtl2api.
#'
#' @return version of qtl2api.
#'
#' @export
version <- function() {
    v <- unlist(utils::packageVersion("qtl2api"))
    paste(v, collapse=".")
}


#' Test if `val` is invalid.
#'
#' @param val value to be tested
#' @return TRUE if invalid, FALSE otherwise
#'
#' @export
invalid <- function(val) {
    return(gtools::invalid(val))
}


#' Test if `val` is valid.
#'
#' @param val value to be tested
#' @return TRUE if valid, FALSE otherwise
#'
#' @export
valid <- function(val) {
    return(!gtools::invalid(val))
}


#' Check value for validity and return it or a default
#'
#' @param value value to check
#' @param default default value to use if value is "invalid"
#'
#' @return `value` if it is valid, `default` otherwise
#'
#' @export
nvl <- function(value, default) {
    if (invalid(value)) {
        return(default)
    }

    value
}


#' Convert value to numeric if possible
#'
#' @param value value to check
#' @param default default value to use if value is "invalid"
#'
#' @return `value` if it is numeric, `default` otherwise
#'
#' @export
nvl_int <- function(value, default) {
    tryCatch(
        {
            n <- as.numeric(value)
            if ((n %% 1) == 0) {
                return(n)
            }
        },
        error = function(cond) {
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    default
}


#' Check dataset to see if the datatype value is "phenotype"
#'
#' @param ds a dataset object
#'
#' @return `TRUE` if the datatype is phenotype, `FALSE` otherwise
#' @export
is_phenotype <- function(ds) {
    if (any("datatype" == names(ds))) {
        if (startsWith(tolower(ds$datatype), "pheno")) {
            return(TRUE)
        }
    } else {
        message("datatype not found in dataset")
    }

    FALSE
}


