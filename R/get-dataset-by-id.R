#' Get the dataset by id (a string).  annot_samples can be used at the top level
#' to share sample annotations among other datasets.
#'
#' @param dataset_id a string, either 'dataset.name' or just 'name'
#'
#' @return the dataset element
#'
#' @export
get_dataset_by_id <- function(dataset_id) {
    if (!is.character(dataset_id)) {
        stop(paste0("dataset_id should be a string, not ", class(dataset_id)))
    }

    dataset <- NULL

    if (exists(dataset_id)) {
        dataset <- get(dataset_id)
    } else {
        expanded <- sprintf("dataset.%s", dataset_id)
        if (exists(expanded)) {
            dataset <- get(expanded)
        }
    }

    if (invalid(dataset)) {
        stop(sprintf("'%s' does not exist", dataset_id))
    }

    return(dataset)
}
