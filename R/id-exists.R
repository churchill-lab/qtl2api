#' Check if id exists and has data in dataset.
#'
#' @param id the id to check
#'
#' @return `TRUE` if id contains data, `FALSE` otherwise
#' @export
id_exists <- function(id) {
    ret <- list()

    datasets <- utils::apropos('^dataset\\.*', ignore.case = TRUE)

    for (d in datasets) {
        ds <- synchronize_dataset(get_dataset_by_id(d))
        all_ids <- ds$annotation_id

        if (id %in% all_ids) {
            ret[[d]] <- list(
                dataset_id        = d,
                dataset_datatype  = ds$datatype,
                dataset_name      = ds$display_name,
                id                = id
            )
        }
    }

    if (length(ret) == 0) {
        ret <- NULL
    }

    ret
}

