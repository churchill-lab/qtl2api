#' Validate the dataset to make sure qtl2api can use it.
#'
#' @param dataset the dataset_id as a string identifier
#' @param extensive TRUE to perform scans on random elements in the data
#'
#' @export
validate_dataset <- function(dataset, extensive = FALSE) {

    # get the dataset
    if (is.character(dataset)) {
        ds_orig <- get_dataset_by_id(dataset)
    } else {
        ds_orig <- dataset
    }

    annots_field <- grep("^display?(\\.|_){1}name$",
                         names(ds_orig),
                         value = TRUE)

    if (length(annots_field) == 0) {
        message("ERROR   : dataset must contain 'display_name'")
    }

    cat("\nSTATUS  : Checking dataset '", ds_orig[[annots_field]], "' ...\n")

    # dataset should be a list
    if (!is.list(ds_orig)) {
        message("ERROR   : dataset should be a list, but found: ", class(ds_orig))
        return()
    }

    if ('datatype' %not in% names(ds_orig)) {
        message("ERROR   : dataset must contain 'datatype'")
    }

    datatype = ds_orig[['datatype']]
    is_mrna <- FALSE
    is_protein <- FALSE
    is_protein_uniprot <- FALSE
    is_phos <- FALSE
    is_pheno <- FALSE
    is_ptm <- FALSE
    is_peptide <- FALSE

    if (tolower(datatype) == 'mrna') {
        is_mrna <- TRUE
    } else if (tolower(datatype) == 'protein') {
        is_protein <- TRUE
    } else if (tolower(datatype) == 'protein_uniprot') {
        is_protein_uniprot <- TRUE
    } else if (tolower(datatype) == 'phos') {
        is_phos <- TRUE
    } else if (tolower(datatype) == 'ptm') {
        is_ptm <- TRUE
    } else if (tolower(datatype) == 'peptide') {
        is_peptide <- TRUE
    } else if (is_phenotype(ds_orig)) {
        is_pheno <- TRUE
    } else {
        message("ERROR   : datatype is invalid: ", datatype)
    }

    validate_annotations(ds_orig)

    ds <- synchronize_dataset(ds_orig)

    if (invalid(ds)) {
        message("ERROR   : dataset not found")
        return()
    }

    validate_covar_info(ds_orig)

    validate_samples(ds_orig)

    validate_lod_peaks(ds_orig)

    validate_data(ds_orig)

    if (extensive) {
        validate_dataset_extensive(ds_orig)
    }
}


#' Perform scans on random data in the dataset designated by dataset_id.
#'
#' @param dataset the dataset as a string identifier
#'
#' @export
validate_dataset_extensive <- function(dataset) {
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

    cat("STATUS  : Checking qtl2api::get_lod_peaks_dataset\n")
    tryCatch(
        {
            # NOTE: dataset instead of ds
            temp <- get_lod_peaks_dataset(dataset)
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

    cat(paste0("STATUS  : Checking qtl2api::get_lod_scan ", id, "\n"))
    tryCatch(
        {
            temp <- get_lod_scan(ds, id)
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
        cat("WARNING : unable to test qtl2api::get_lod_scan with intcovar\n")
    } else {
        cat(paste0("STATUS  : Checking qtl2api::get_lod_scan: ", id, " intcovar: ", intcovar, "\n"))
        tryCatch(
            {
                temp <- get_lod_scan(ds, id, intcovar)
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
                temp <- get_lod_scan_by_sample(ds, id, chrom, intcovar)
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

    cat(paste0("STATUS  : Checking qtl2api::get_founder_coefficients: ", id, " chrom: ", chrom, "\n"))
    tryCatch(
        {
            temp <- get_founder_coefficients(ds, id, chrom)
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
        cat("WARNING : unable to test qtl2api::get_founder_coefficients with intcovar\n")
    } else {
        cat(paste0("STATUS  : Checking qtl2api::get_founder_coefficients: ", id, " chrom: ", chrom, ", intcovar: ", intcovar, "\n"))
        tryCatch(
            {
                temp <- get_founder_coefficients(ds, id, chrom, intcovar = intcovar)
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

    cat(paste0("STATUS  : Checking qtl2api::get_correlation ", id, "\n"))
    tryCatch(
        {
            temp <- get_correlation(ds, id)
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
        cat("WARNING : unable to test qtl2api::get_correlation with intcovar\n")
    } else {
        cat(paste0("STATUS  : Checking qtl2api::get_correlation ", id, " intcovar: ", intcovar, "\n"))
        tryCatch(
            {
                temp <- get_correlation(ds, id, intcovar = intcovar)
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
    cat(paste0("STATUS  : Checking qtl2api::get_correlation_plot_data ", id, " id_correlate: ", id_correlate, "\n"))
    tryCatch(
        {
            temp <- get_correlation_plot_data(ds, id, ds, id_correlate)
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
        cat("WARNING : unable to test qtl2api::get_correlation_plot_data with intcovar\n")
    } else {
        cat(paste0("STATUS  : Checking qtl2api::get_correlation_plot_data ", id, " id_correlate: ", id_correlate, " intcovar: ", intcovar, "\n"))
        tryCatch(
            {
                temp <- get_correlation_plot_data(ds, id, ds, id_correlate, intcovar = intcovar)
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
        cat(paste0("STATUS  : Checking qtl2api::get_mediation ", id, " marker: ", marker, "\n"))
        tryCatch(
            {
                temp <- get_mediation(ds, id, marker)
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
}
