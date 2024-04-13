#' Get the samples from the dataset.
#'
#' @param dataset the dataset object
#'
#' @return the samples tibble
#'
#' @export
get_samples <- function(dataset) {
    # allow some play in the name sample, samples
    samples_name <- grep("^samples?$", names(dataset), value = TRUE)

    if (length(samples_name) == 0) {
        stop("Unable to find any samples tibble in dataset")
    }

    samples <-
        dataset[[samples_name]]

    # again allow some play in the column name, even allow for mouse
    samples_id_field <- grep(
        "^mouse$|^sample$|^mouse(\\.|_){0,1}?id$|^samples?(\\.|_){0,1}?id$",
        colnames(samples),
        value = TRUE,
        ignore.case = TRUE
    )

    # make sure we pass back a tibble with sample_id
    if (samples_id_field != "sample_id") {
        samples$sample_id <- samples[[samples_id_field]]
    }

    # make sure samples_id is the first column (easier to view)
    samples <- samples %>% dplyr::relocate(sample_id = .data$sample_id)

    return (samples)
}
