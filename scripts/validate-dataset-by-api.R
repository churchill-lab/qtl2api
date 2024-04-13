#' Perform scans on random data in the dataset designated by dataset_id.
#'
#' @param dataset the dataset as a string identifier
validate_dataset_by_api <- function(dataset) {
    cat("STATUS  : Checking qtl2api operations...\n")

    # get the dataset
    if (is.character(dataset)) {
        ds_orig <- get_dataset_by_id(dataset)
    } else {
        ds_orig <- dataset
    }

    ds <- synchronize_dataset(ds_orig)

    datatype = ds[['datatype']]
    id <- get_random_id(ds)

    intcovar <- NULL

    tryCatch(
        {
            covar_temp <-
                ds$covar_info %>%
                dplyr::filter(.data$primary == TRUE)

            intcovar <- covar_temp$sample_column
            # TODO: Loop through all covars
            intcovar <- intcovar[1]
        },
        error = function(cond) {
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    chrom <- sample(names(K), 1)

    markers_temp <- markers %>%
        dplyr::sample_n(1)

    marker <- markers_temp$marker.id



    cat("STATUS  : Checking qtl2api::get_lod_peaks\n")
    tryCatch(
        {
            # NOTE: dataset instead of ds
            temp <- get_lod_peaks(dataset)
        },
        error = function(cond) {
            message("ERROR   : ", cond$message)
        },
        warning = function(cond) {
            message("ERROR   : ", cond$message)
        },
        finally = {
        }
    )

    cat(paste0("STATUS  : Checking qtl2api::calc_lod_scores ", id, "\n"))
    tryCatch(
        {
            temp <- calc_lod_scores(ds, id)
        },
        error = function(cond) {
            message("ERROR   : ", cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    if (invalid(intcovar)) {
        cat("WARNING : unable to test qtl2api::calc_lod_scores with intcovar\n")
    } else {
        cat(paste0("STATUS  : Checking qtl2api::calc_lod_scores: ", id, " intcovar: ", intcovar, "\n"))
        tryCatch(
            {
                temp <- calc_lod_scores(ds, id, intcovar)
            },
            error = function(cond) {
                message("ERROR   : ", cond$message)
            },
            warning = function(cond) {
            },
            finally = {
            }
        )
    }

    if (invalid(intcovar)) {
        cat("WARNING : unable to test qtl2api::get_lod_scan_by_sample with intcovar\n")
    } else {
        cat(paste0("STATUS  : Checking qtl2api::get_lod_scan_by_sample: ", id, " chrom: ", chrom, " intcovar: ", intcovar, "\n"))
        tryCatch(
            {
                temp <- calc_lod_scores_by_covar(ds, id, chrom, intcovar)
            },
            error = function(cond) {
                message("ERROR   : ", cond$message)
            },
            warning = function(cond) {
            },
            finally = {
            }
        )
    }

    cat(paste0("STATUS  : Checking qtl2api::get_expression: ", id, "\n"))
    tryCatch(
        {
            temp <- get_expression(ds, id)
        },
        error = function(cond) {
            message("ERROR   : ", cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    cat(paste0("STATUS  : Checking qtl2api::calc_founder_coefficients: ", id, " chrom: ", chrom, "\n"))
    tryCatch(
        {
            temp <- calc_founder_coefficients(ds, id, chrom)
        },
        error = function(cond) {
            message("ERROR   : ", cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    if (invalid(intcovar)) {
        cat("WARNING : unable to test qtl2api::calc_founder_coefficients with intcovar\n")
    } else {
        cat(paste0("STATUS  : Checking qtl2api::calc_founder_coefficients: ", id, " chrom: ", chrom, ", intcovar: ", intcovar, "\n"))
        tryCatch(
            {
                temp <- calc_founder_coefficients(ds, id, chrom, intcovar = intcovar)
            },
            error = function(cond) {
                message("ERROR   : ", cond$message)
            },
            warning = function(cond) {
            },
            finally = {
            }
        )
    }

    cat(paste0("STATUS  : Checking qtl2api::calc_correlation ", id, "\n"))
    tryCatch(
        {
            temp <- calc_correlation(ds, id)
        },
        error = function(cond) {
            message("ERROR   : ", cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    if (invalid(intcovar)) {
        cat("WARNING : unable to test qtl2api::calc_correlation with intcovar\n")
    } else {
        cat(paste0("STATUS  : Checking qtl2api::calc_correlation ", id, " intcovar: ", intcovar, "\n"))
        tryCatch(
            {
                temp <- calc_correlation(ds, id, intcovar = intcovar)
            },
            error = function(cond) {
                message("ERROR   : ", cond$message)
            },
            warning = function(cond) {
            },
            finally = {
            }
        )
    }

    id_correlate <- get_random_id(ds)
    cat(paste0("STATUS  : Checking qtl2api::calc_correlation_plot ", id, " id_correlate: ", id_correlate, "\n"))
    tryCatch(
        {
            temp <- calc_correlation_plot(ds, id, ds, id_correlate)
        },
        error = function(cond) {
            message("ERROR   : ", cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    if (invalid(intcovar)) {
        cat("WARNING : unable to test qtl2api::calc_correlation_plot with intcovar\n")
    } else {
        cat(paste0("STATUS  : Checking qtl2api::calc_correlation_plot ", id, " id_correlate: ", id_correlate, " intcovar: ", intcovar, "\n"))
        tryCatch(
            {
                temp <- calc_correlation_plot(ds, id, ds, id_correlate, intcovar = intcovar)
            },
            error = function(cond) {
                message("ERROR   : ", cond$message)
            },
            warning = function(cond) {
            },
            finally = {
            }
        )
    }

    if (!is_phenotype(ds)) {
        cat(paste0("STATUS  : Checking qtl2api::calc_mediation ", id, " marker: ", marker, "\n"))
        tryCatch(
            {
                #temp <-
                    calc_mediation(ds, id, marker)
            },
            error = function(cond) {
                message("ERROR   : ", cond$message)
            },
            warning = function(cond) {
                message("WARNING  : ", cond$message)
            },
            finally = {
            }
        )
    }
}
