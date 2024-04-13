#' Validate the annot_samples in the dataset to make sure qtl2api can use it.
#'
#' @param dataset the dataset as a string identifier
#'
#' @export
validate_samples <- function(dataset) {
    cat("STATUS  : Checking samples\n")

    # get the dataset
    if (is.character(dataset)) {
        ds_orig <- get_dataset_by_id(dataset)
    } else {
        ds_orig <- dataset
    }

    ds <- synchronize_dataset(ds_orig)

    samples_name <- grep("^samples?$", names(ds_orig), value = TRUE)

    if (length(samples_name) == 0) {
        message('ERROR   : samples not found')
        return()
    }

    if (!tibble::is_tibble(ds_orig[[samples_name]])) {
        message(paste0("ERROR   : samples should be a tibble, but found ", class(ds[[samples_name]])))
    }

    if (any(duplicated(ds[[samples_name]]["sample_id"]))) {
        message("ERROR   : there are duplicated annotations in samples")
    }

    num_annots_orig <- NROW(ds_orig[[samples_name]])
    num_annots_synch <- NROW(ds$samples)

    if (num_annots_orig != num_annots_synch) {
        cat("WARNING : synchronizing samples, # samples changed from", num_annots_orig, "to", num_annots_synch, "\n")
    }

    # let's see which columns are not FACTORS, but should be
    for (x in ds$covar_info$sample_column) {
        n <- length(unique(ds$samples[[x]]))
        if (is.factor(ds$samples[[x]])) {
            if (n == 1) {
                cat("WARNING : ", x, "is listed as a covar and a factor class, but only 1 unique value (remove from covar info?)\n")
            }
        } else {
            cat("WARNING : ", x, "is listed as a covar, not factor class,", n, "unique values found (change to factor?)\n")
        }
    }
}
