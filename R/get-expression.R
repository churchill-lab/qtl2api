#' Get the expression data.
#'
#' @param dataset the dataset object
#' @param id the unique id in the dataset
#'
#' @return a `list` with elements of tibble and list of the datatypes.
#'
#' @importFrom rlang .data
#' @export
get_expression <- function(dataset, id) {
    # make sure annotations, data, and samples are synchronized
    ds <- synchronize_data(dataset)

    # check if id exists
    idx <- which(colnames(ds$data) == id)

    if (gtools::invalid(idx)) {
        stop(sprintf("Cannot find data for '%s' in dataset", id))
    }

    if (gtools::invalid(ds$covar_info)) {
        datatypes <- NULL
    } else {
        datatypes <- list()
        for (i in ds$covar_info$sample_column) {

            stopifnot(!is.null(ds$annot_samples[[i]]))

            if (is.factor(ds$annot_samples[[i]])) {
                datatypes[[i]] <-
                    gtools::mixedsort(levels(ds$annot_samples[[i]]))
            } else {
                datatypes[[i]] <-
                    gtools::mixedsort(unique(ds$annot_samples[[i]]))
            }
        }
    }

    # get the sample id field
    sample_id_field <- get_sample_id_field(ds)

    # only pass back the columns we need
    samples <- ds$annot_samples %>%
        dplyr::select(
            sample_id = sample_id_field,
            dplyr::all_of(names(datatypes))
        )

    # bind the data
    expression_temp <- tibble::tibble(
        sample_id  = rownames(ds$data),
        expression = ds$data[, idx]
    )

    output <- samples %>%
        dplyr::inner_join(
            expression_temp,
            by = "sample_id"
        )

    list(
        data = output,
        datatypes = datatypes
    )
}

