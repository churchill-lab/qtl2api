#' Get the data from a data from the dataset.
#'
#' The data element in the dataset list should either be a matrix or a named
#' list with each eleemnt being a matrix.
#'
#' If `data_name` is not specified, the order that will be returned is:
#'     data (if a matrix)
#'     rz, norm, raw, log, transformed (if a named list).
#'
#' @param ds a dataset object (synchronized or not)
#' @param data_name A string containing which data element from the dataset's
#'     data element
#'
#' @return the data element
#' @export
get_data <- function(ds, data_name = NULL) {
    ret <- NULL

    # check if data is a matrix, but something else was requested
    if (is.matrix(ds$data) && (!is.null(data_name))) {
        stop(sprintf("Specified data '%s' not found in dataset", data_name))
    }

    if (valid(data_name)) {
        if (is.matrix(ds$data)) {
            stop(sprintf("Specified data '%s' not found in dataset", data_name))
        }

        # return the requested data element in the named list
        ret <- ds$data[[data_name]]
    } else {
        # Order of return is data (if matrix), than:
        # rz, norm, raw, log, transformed
        if (is.matrix(ds$data)) {
            ret <- ds$data
        } else {
            if (!is.null(ds$data$rz)) {
                ret <- ds$data$rz
            } else if (!is.null(ds$data$norm)) {
                ret <- ds$data$norm
            } else if (!is.null(ds$data$raw)) {
                ret <- ds$data$raw
            } else if (!is.null(ds$data$log)) {
                ret <- ds$data$log
            } else if (!is.null(ds$data$transformed)) {
                ret <- ds$data$transformed
            }
        }
    }

    if (invalid(ret)) {
        stop("Unable to find data in dataset")
    }

    ret
}
