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

    annots_field <- grep("^annots?(\\.|_){1}samples?$",
                         names(ds_orig),
                         value = TRUE)

    if (length(annots_field) == 0) {
        message('ERROR   : annot_samples not found')
        return()
    }

    if (!tibble::is_tibble(ds_orig[[annots_field]])) {
        message(paste0("ERROR   : annot_samples should be a tibble, but found ", class(ds[[annots_field]])))
    }

    sample_id_field <- ds$sample_id_field
    if (invalid(sample_id_field)) {
        message('ERROR   : unable to determine sample id field')
    }

    if (any(duplicated(ds[[annots_field]][sample_id_field]))) {
        message("ERROR   : there are duplicated annotations in annot_samples")
    }

    num_annots_orig <- NROW(ds_orig[[annots_field]])
    num_annots_synch <- NROW(ds$annot_samples)

    if (num_annots_orig != num_annots_synch) {
        cat("WARNING : synchronizing samples, # samples changed from", num_annots_orig, "to", num_annots_synch, "\n")
    }

    # let's see which columns are not FACTORS, but should be
    for (x in ds$covar_info$sample_column) {
        n <- length(unique(ds$annot_samples[[x]]))
        if (is.factor(ds$annot_samples[[x]])) {
            if (n == 1) {
                cat("WARNING : ", x, "is listed as a covar and a factor class, but only 1 unique value (remove from covar info?)\n")
            }
        } else {
            cat("WARNING : ", x, "is listed as a covar, not factor class,", n, "unique values found (change to factor?)\n")
        }
    }
}
