#' Get the expression data.
#'
#' @param ds the dataset object
#' @param id the unique id in the dataset
#'
#' @return a `list` with elements of tibble and list of the datatypes.
#'
#' @importFrom rlang .data
#' @export
get_expression <- function(ds, id) {
    # get the data
    data <- get_data(ds)

    # check if id exists
    idx <- which(colnames(data) == id)

    if (gtools::invalid(idx)) {
        stop(sprintf("Cannot find data for '%s' in dataset", id))
    }

    if (gtools::invalid(ds$covar.info)) {
        datatypes <- NULL
    } else {
        datatypes <- list()
        covar_info <- ds$covar.info %>% janitor::clean_names()

        for (i in covar_info$sample_column) {
            stopifnot(!is.null(ds$annot.samples[[i]]))
            if (is.factor(ds$annot.samples[[i]])) {
                datatypes[[i]] <-
                    gtools::mixedsort(levels(ds$annot.samples[[i]]))
            } else {
                datatypes[[i]] <-
                    gtools::mixedsort(unique(ds$annot.samples[[i]]))
            }
        }
    }

    # get the sample id field
    sample_id_field <- get_sample_id_field(ds)

    # only pass back the columns we need
    samples <- ds$annot.samples %>%
        dplyr::select(
            sample_id = sample_id_field,
            dplyr::all_of(names(datatypes))
        )

    # bind the data
    expression_temp <- tibble::tibble(
        sample_id  = rownames(data),
        expression = data[, idx]
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
