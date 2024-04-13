#' Validate the covar_info in the dataset to make sure qtl2api can use it.
#'
#' @param dataset the dataset_id as a string identifier
#'
#' @export
validate_covar_info <- function(dataset) {
    cat("STATUS  : Checking covar_info\n")

    # get the dataset
    if (is.character(dataset)) {
        ds_orig <- get_dataset_by_id(dataset)
    } else {
        ds_orig <- dataset
    }

    ds <- synchronize_dataset(ds_orig)

    annots_field <- grep("^covar(\\.|_){1}info$",
                         names(ds_orig),
                         value = TRUE)

    if (length(annots_field) == 0) {
        message("WARNING : covar_info not found")
        return()
    }

    if (is.null(ds_orig[[annots_field]])) {
        message("WARNING : covar_info is NULL")
        return()
    }

    if (!tibble::is_tibble(ds_orig[[annots_field]])) {
        message("covar_info should be a tibble, but found '", class(ds_orig[[annots_field]]), "'")
    }

    column_field <- grep(
        "^sample(\\.|_){1}column$",
        colnames(ds_orig[[annots_field]]),
        value = TRUE,
        ignore.case = TRUE
    )

    if (length(column_field) == 0) {
        message("ERROR   : sample_column not found in covar_info")
    }

    column_field <- grep(
        "^display(\\.|_){1}name$",
        colnames(ds_orig[[annots_field]]),
        value = TRUE,
        ignore.case = TRUE
    )

    if (length(column_field) == 0) {
        message("ERROR   : display_name not found in covar_info")
    }

    column_field <- grep(
        "^lod(\\.|_){1}peaks$",
        colnames(ds_orig[[annots_field]]),
        value = TRUE,
        ignore.case = TRUE
    )

    if (length(column_field) == 0) {
        message("ERROR   : lod_peaks not found in covar_info")
    }

    if ('interactive' %not in% names(ds_orig[[annots_field]])) {
        message("ERROR   : interactive not found in covar_info")
    } else if ('primary' %not in% names(ds_orig[[annots_field]])) {
        message("ERROR   : primary not found in covar_info")
    }

    if(!is.logical(ds_orig[[annots_field]]$primary)) {
        message('ERROR   : covar_info$primary should be logical, not ', class(ds_orig[[annots_field]]$primary))
    }

    if(!is.logical(ds_orig[[annots_field]]$interactive)) {
        message('ERROR   : covar_info$interactive should be logical, not ', class(ds_orig[[annots_field]]$interactive))
    }

    for(i in 1:nrow(ds$covar_info)) {
        row <- ds$covar_info[i, ]

        if (row$sample_column %not in% colnames(ds$samples)) {
            message("ERROR   : covar_info$sample_column ('", row$sample_column, "') is not a column name in samples")
        }

        if (invalid(row$display_name)) {
            message("ERROR   : covar_info$display_name needs to have a value")
        }

        if (class(row$interactive) == 'logical') {
            if (row$interactive) {
                if (invalid(row$lod_peaks)) {
                    message("ERROR   : covar_info$interactive is TRUE, but covar_info$lod_peaks is NA")
                } else {
                    # check for existence of lod_peaks
                    annots_field_peaks <- grep("^lod(\\.|_){1}peaks?$",
                                               names(ds_orig),
                                               value = TRUE)

                    if ((length(annots_field_peaks) != 0) && (!is.null(ds_orig[[annots_field_peaks]]))) {
                        if (row$lod_peaks %not in% names(ds_orig[[annots_field_peaks]])) {
                            message("ERROR   : covar_info$interactive is TRUE, but covar_info$lod_peaks ('", row$lod_peaks, "') is not in lod_peaks")
                        }
                    }
                }
            } else {
                if (valid(row$lod_peaks)) {
                    message("ERROR   : covar_info$interactive is FALSE, but covar_info$lod_peaks ('", row$lod_peaks, "') is set")
                }
            }
        } else {
            message("ERROR   : covar_info$interactive is INVALID, please set correctly")
        }
    }

    if (class(ds$covar_info$primary) == "logical") {
        if (!any(ds$covar_info$primary)) {
            message("ERROR   : covar_info$primary needs 1 value set to TRUE")
        }
    } else {
        message("ERROR   : covar_info$primary class is not logical (TRUE or FALSE), please set correctly")
    }

    cat("STATUS  : Checking covar_matrix\n")

    if (is_phenotype(ds)) {
        phenos <- ds$annotations %>% dplyr::filter(.data$is_pheno == TRUE)
        phenos <- phenos$data_name

        i <- 0
        for (x in phenos) {
            i <- i + 1
            tryCatch(
                {
                    if (i %% 500 == 0) {
                        cat("STATUS  : ... completed", i, "checks out of", length(phenos), "\n")
                    }
                    temp <- get_covar_matrix(ds, x)
                },
                error = function(cond) {
                    message('ERROR   : Check phenotype ', x)
                    message('ERROR   : Unable to generate covar_matrix, check covar_info')
                    message(cond$message)
                },
                warning = function(cond) {
                },
                finally = {
                }
            )
        }

    } else {
        tryCatch(
            {
                temp <- get_covar_matrix(ds)
            },
            error = function(cond) {
                message('ERROR   : Unable to generate covar_matrix, check covar_info')
                message(cond$message)
            },
            warning = function(cond) {
            },
            finally = {
            }
        )
    }
}
