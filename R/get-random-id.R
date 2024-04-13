#' Get a random identifier from the dataset.
#'
#' @param dataset a dataset object
#'
#' @return a random identifier
#' @export
get_random_id <- function(dataset) {
    if (invalid(dataset)) {
        stop("invalid dataset object")
    } else if (!is.list(dataset)) {
        stop(paste0("dataset should be a list, but it's a ", class(dataset)))
    }

    # make sure the data has been synchronized
    if (invalid(dataset$is_synchronized)) {
        dataset <- synchronize_dataset(dataset)
    }

    return(sample(dataset$annotations$annotation_id, 1))
}
