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

    ds <- synchronize_dataset(ds_orig)

    if (gtools::invalid(ds)) {
        message("ERROR   : dataset not found")
        return()
    }

    # dataset should be a list
    if (!is.list(ds_orig)) {
        message("ERROR   : dataset should be a list, but found: ", class(ds_orig))
        return()
    }

    if ('datatype' %not in% names(ds_orig)) {
        message("ERROR   : dataset must contain 'datatype'")
    }

    datatype = ds[['datatype']]
    is_mrna <- FALSE
    is_protein <- FALSE
    is_phos <- FALSE
    is_pheno <- FALSE

    if (tolower(datatype) == 'mrna') {
        is_mrna <- TRUE
    } else if (tolower(datatype) == 'protein') {
        is_protein <- TRUE
    } else if (tolower(datatype) == 'phos') {
        is_phos <- TRUE
    } else if (is_phenotype(ds)) {
        is_pheno <- TRUE
    } else {
        message("ERROR   : datatype is invalid: ", datatype)
    }

    validate_annotations(ds_orig)

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

    if (gtools::invalid(intcovar)) {
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

    if (gtools::invalid(intcovar)) {
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

    if (gtools::invalid(intcovar)) {
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

    if (gtools::invalid(intcovar)) {
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

    if (gtools::invalid(intcovar)) {
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


#' Validate the environment for defined objects.
#'
#' @param extensive TRUE to perform scans on random elements in the data
#'
#' @export
validate_environment <- function(extensive = FALSE) {
    cat("\n")
    cat("qtl2api tries to use gene_id/gene.id, protein_id/protein.id,\n")
    cat("data_name/data.name along as trying to figure out annotations for\n")
    cat("mouse_id/mouse.id/sample_id/sample.id\n")
    cat("\n")
    cat("Internally, qtl2api uses underscore '_' for variable names and\n")
    cat("column names\n")
    cat("\n")
    cat("If an error is detected, this method will use the '_' version\n")
    cat("\n")


    # grab the datasets in the environment
    datasets <- grep('^dataset*', utils::apropos('dataset\\.'), value=TRUE)

    if (gtools::invalid(datasets)) {
        message("ERROR   : No datasets found!")
        return()
    }

    ensembl_version <-
        utils::apropos("^ensembl(\\.|_){1}version$", ignore.case = TRUE)

    if (length(ensembl_version) == 0) {
        message("ERROR   : ensembl_version does not exist")
    }

    # Check genoprobs and K.
    if(length(genoprobs) != length(K)) {
        message("ERROR   : genoprobs (", length(genoprobs), ") and K (", length(K), ") do not have the same length.")
    } else {
        if(any(names(genoprobs) != names(K))) {
            message("ERROR   : names of genoprobs and K do not match.")
        }

        rownames.eq <- mapply(function(x, y) { all(rownames(x) == rownames(y)) }, genoprobs, K)
        if(any(rownames.eq == FALSE)) {
            message("ERROR   : sample IDs do not match between genoprobs and K.")
        }
    }

    # Check Marker IDs for genoprobs and map
    if(length(genoprobs) != length(map)) {
        message("ERROR   : genoprobs (", length(genoprobs), ") and map (", length(map), ") do not have the same length.")
    } else {
        rownames.eq <- mapply(function(x, y) { all(dimnames(x)[[3]] == names(y)) }, genoprobs, map)
        if(any(rownames.eq == FALSE)) {
            message("ERROR   : marker names do not match between genoprobs and map")
        }
    }

    # Check dimensions of markers and map.
    map.length = sum(sapply(map, length))
    if(map.length != nrow(markers)) {
        message("ERROR   : number of rows in markers (", nrow(markers), ") does not equal the number of markers in map (", map.length, ")")
    }

    for (ds in datasets) {
        validate_dataset(ds, extensive)
    }
}


#' Validate the annotations in the dataset to make sure qtl2api can use it.
#'
#' @param dataset the dataset as a string identifier
#'
#' @export
validate_annotations <- function(dataset) {
    cat("STATUS  : Checking annotations\n")

    # get the dataset
    if (is.character(dataset)) {
        ds_orig <- get_dataset_by_id(dataset)
    } else {
        ds_orig <- dataset
    }

    ds <- synchronize_dataset(ds_orig)

    annots <- NULL
    annots_orig <- NULL

    if (tolower(ds$datatype) == "mrna") {
        annots_field <- grep("^annots?(\\.|_){1}mrnas?$",
                             names(ds_orig),
                             value = TRUE)

        if (length(annots_field) == 0) {
            message("ERROR   : annot_mrna not found in dataset.")
            return()
        }

        if (!tibble::is_tibble(ds_orig[[annots_field]])) {
            message("ERROR   : annot_mrna should be a tibble, but found: ", class(ds[[annots_field]]))
        }

        annots <- ds$annot_mrna
        annots_orig <- ds_orig[[annots_field]]

        if (any(duplicated(annots$gene_id))) {
            message("ERROR   : There are duplicated gene identifiers in annot_mrna")
        }
    } else if (tolower(ds$datatype) == "protein") {
        annots_field <- grep("^annots?(\\.|_){1}proteins?$",
                             names(ds_orig),
                             value = TRUE)

        if (length(annots_field) == 0) {
            message("ERROR   : annot_protein not found in dataset.")
            return()
        }

        if (!tibble::is_tibble(ds_orig[[annots_field]])) {
            message("ERROR   : annot_protein should be a tibble, but found: ", class(ds[[annots_field]]))
        }

        annots <- ds$annot_protein
        annots_orig <- ds_orig[[annots_field]]

        if (any(duplicated(annots$protein_id))) {
            message("ERROR   : There are duplicated protein identifiers in annot_protein")
        }
    } else if (tolower(ds$datatype) == "phos") {
        annots_field <- grep("^annots?(\\.|_){1}phos$",
                             names(ds_orig),
                             value = TRUE)

        if (length(annots_field) == 0) {
            message("ERROR   : annot_phos not found in dataset.")
            return()
        }

        if (!tibble::is_tibble(ds_orig[[annots_field]])) {
            message("ERROR   : annot_phos should be a tibble, but found: ", class(ds[[annots_field]]))
        }

        annots <- ds$annot_phos
        annots_orig <- ds_orig[[annots_field]]

        if (any(duplicated(annots$phos_id))) {
            message("ERROR   : There are duplicated phos identifiers in annot_phos")
        }
    } else if (is_phenotype(ds)) {
        annots_field <- grep("^annots?(\\.|_){1}pheno(type)?s?$",
                             names(dataset),
                             value = TRUE)

        if (length(annots_field) == 0) {
            message("ERROR   : annot_phenotype not found in dataset.")
            return()
        }

        if (!tibble::is_tibble(ds_orig[[annots_field]])) {
            message("ERROR   : annot_phenotype should be a tibble, but found: ", class(ds[[annots_field]]))
        }

        annots <- ds$annot_phenotype
        annots_orig <- ds_orig[[annots_field]]

        if (any(duplicated(annots$data_name))) {
            message("ERROR   : There are duplicated data_name annotations in annot_phenotype")
        }
    }

    if (is_phenotype(ds)) {
        expected_names <- c('data_name', 'short_name', 'description', 'is_id',
                            'category', 'is_numeric', 'is_date', 'is_factor',
                            'factor_levels', 'is_pheno', 'omit', 'use_covar')

        for (n in expected_names) {
            if (n %not in% names(annots)) {
                message("ERROR   :", n, "not found in annot_phenotype")
            }
        }

        expected_logical <- c('is_id', 'is_pheno', 'is_numeric', 'is_date',
                              'is_factor', 'omit')

        for (n in expected_logical) {
            if(!is.logical(annots[[n]])) {
                message("ERROR   : annot_phenotype$", n, " should be logical, not ", class(annots[[n]]))
            }
        }

        column_field <- grep(
            "^is(\\.|_){1}id$",
            colnames(annots_orig),
            value = TRUE,
            ignore.case = TRUE
        )

        column_field_data_name <- grep(
            "^data(\\.|_){1}name$",
            colnames(annots_orig),
            value = TRUE,
            ignore.case = TRUE
        )

        if(is.logical(annots_orig[[column_field]])) {
            # one and only 1 ID
            theID <- annots_orig[which(annots_orig[[column_field]] == TRUE),][[column_field_data_name]]

            if(length(theID) != 1) {
                message("ERROR   : annot_phenotype$is_id should have 1 and only 1 row set to TRUE")
            }
        }
    } else {
        annot_name <- NULL

        if (tolower(ds$datatype) == "mrna") {
            annot_name <- "annot_mrna"
        } else if (tolower(ds$datatype) == "protein") {
            annot_name <- 'annot_protein'

            if ('protein_id' %not in% names(annots)) {
                message("ERROR   : protein_id not found in annot_protein")
            }
        } else if (tolower(ds$datatype) == "phos") {
            annot_name <- 'annot_phos'

            if ('protein_id' %not in% names(annots)) {
                message("ERROR   : protein_id not found in annot_phos")
            }

            annot_name <- 'annot_phos'

            if ('phos_id' %not in% names(annots)) {
                message("ERROR   : phos_id not found in annot_phos")
            }
        }

        expected_names <-  c('gene_id', 'symbol', 'chr', 'start', 'end')

        for (n in expected_names) {
            if (n %not in% names(annots)) {
                message("ERROR   : ", n, " not found in ", annot_name)
            }
        }

        #if (any(annots$start > 10000.0)) {
        #    message(annot_name, '$start should be in Mbp not bp')
        #} else if (any(annots$end > 10000.0)) {
        #    message(annot_name, '$end should be in Mbp not bp')
        #}
    }

    num_annots_orig <- NROW(annots_orig)
    num_annots_synch <- NROW(annots)

    if (num_annots_orig != num_annots_synch) {
        cat("WARNING : synchronizing annotations, # samples changed from", num_annots_orig, "to", num_annots_synch, "\n")
    }
}


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
    if (gtools::invalid(sample_id_field)) {
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

        if (row$sample_column %not in% colnames(ds$annot_samples)) {
            message("ERROR   : covar_info$sample_column ('", row$sample_column, "') is not a column name in annot_samples")
        }

        if (gtools::invalid(row$display_name)) {
            message("ERROR   : covar_info$display_name needs to have a value")
        }

        if (class(row$interactive) == 'logical') {
            if (row$interactive) {
                if (gtools::invalid(row$lod_peaks)) {
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
                if (!gtools::invalid(row$lod_peaks)) {
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
        phenos <- ds$annot_phenotype %>% dplyr::filter(.data$is_pheno == TRUE)
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


#' Validate the lod.peaks to make sure qtl2api can use it.
#'
#' @param dataset the dataset as a string identifier
#'
#' @export
validate_lod_peaks <- function(dataset) {
    cat("STATUS  : Checking lod_peaks\n")

    # get the dataset
    if (is.character(dataset)) {
        ds_orig <- get_dataset_by_id(dataset)
    } else {
        ds_orig <- dataset
    }

    ds <- synchronize_dataset(ds_orig)

    annots_field <- grep("^lod(\\.|_){1}peaks$",
                         names(ds_orig),
                         value = TRUE)

    if (length(annots_field) == 0) {
        message("WARNING : lod_peaks not found")
        return()
    }

    if (!is.list(ds_orig[[annots_field]])) {
        message("ERROR   : lod_peaks should be a list, but found ", class(ds_orig[[annots_field]]))
    }

    if ('additive' %not in% names(ds_orig[[annots_field]])) {
        message("ERROR   : additive should be an element in lod_peaks")
    }

    # get the rest
    for (i in seq(nrow(ds$covar_info))) {
        cov_inf <- ds$covar_info[i, ]

        # only look at interactive peaks
        if (cov_inf$interactive) {
            annots_field_peaks <- grep("^lod(\\.|_){1}peaks?$",
                                       names(ds_orig),
                                       value = TRUE)

            peaks <- ds_orig[[annots_field_peaks]][[cov_inf$lod_peaks]]

            if (gtools::invalid(peaks)) {
                message("ERROR   : Unable to find lod_peaks '", cov_inf$lod_peaks, "'")
            } else {

                peaks <- peaks %>%
                    janitor::clean_names()

                if (tolower(ds$datatype) == "protein") {
                    if (length(setdiff(peaks$protein_id, ds$annot_protein$protein_id))) {
                        cat("WARNING : not all lod_peaks$protein_id are in annot_protein$protein_id\n")
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        cat("WARNING : not all lod_peaks$marker_id are in markers\n")
                    }
                } else if (tolower(ds$datatype) == "mrna") {
                    if (length(setdiff(peaks$gene_id, ds$annot_mrna$gene_id))) {
                        cat("WARNING : not all lod_peaks$gene_id are in annot_mrna$gene_id\n")
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        cat("WARNING : not all lod_peaks$marker_id are in markers\n")
                    }
                } else {
                    if (length(setdiff(peaks$data_name, ds$annot_phenotype$data_name))) {
                        cat("WARNING : not all lod_peaks['", cov_inf$sample_column, "']$data_name are in annot_phenotype$data_name\n")
                    }
                }
            }
        }
    }
}

#' Validate the data in the dataset to make sure qtl2api can use it.
#'
#' @param dataset the dataset_id as a string identifier
#'
#' @export
validate_data <- function(dataset) {
    cat("STATUS  : Checking data\n")

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

    num_annot_samples <- nrow(ds_orig[[annots_field]])

    if (tolower(ds$datatype) == "mrna") {
        annot_ids <- ds$annot_mrna$gene_id
    } else if (tolower(ds$datatype) == "protein") {
        annot_ids <- ds$annot_protein$protein_id
    } else if (tolower(ds$datatype) == "phos") {
        annot_ids <- ds$annot_phos$phos_id
    } else {
        annots_temp <-
            ds$annot_phenotype %>%
            dplyr::filter(
                .data$omit == FALSE,
                .data$is_pheno == TRUE
            )
        annot_ids <- annots_temp$data_name
    }

    if ("data" %not in% names(ds)) {
        message("ERROR   : data not found")
        return()
    }

    if (is.matrix(ds$data)) {
        # check if the data is numeric
        if (!is.numeric(ds$data)) {
            message("ERROR   : data matrix is not numeric")
            return()
        }

        if (num_annot_samples != nrow(ds$data)) {
            cat("WARNING : number of samples (", num_annot_samples, ") != ",
                "number of data rows (", nrow(ds$data), ")\n")
        }

        if (NCOL(ds$data) != NROW(annot_ids)) {
            cat('WARNING : number of annotations (', NROW(annot_ids), ') != ',
                'number of data cols (', NCOL(ds$data), ")\n")

            x <- setdiff(annot_ids, colnames(ds$data))
            y <- setdiff(colnames(ds$data), annot_ids)

            if (length(x) > 0) {
                #cat("WARNING : annotations with no data:", paste(x, sep = ",", collapse = ","), "\n")
                cat("WARNING : # annotations with no data:", length(x), "\n")
            }

            if (length(y) > 0) {
                #cat("WARNING : data with no annotations:", paste(y, sep = ",", collapse = ","), "\n")
                cat("WARNING : # data with no annotations:", length(y), "\n")
            }
        }

        #perc_missing = (sum(is.na(ds_orig$data)) * 100) / prod(dim(ds_orig$data))
        #cat("STATUS  : Percentage of missing original data: ", perc_missing, "\n")

        perc_missing = (sum(is.na(ds$data)) * 100) / prod(dim(ds$data))
        if (perc_missing >= 10) {
            cat("WARNING : Percentage of missing synchronized data: ", perc_missing, "\n")
        } else {
            cat("STATUS  : Percentage of missing synchronized data: ", perc_missing, "\n")
        }

    } else if (is.list(ds$data)) {
        # TODO: this won't happen on a synchronized dataset
        # no list, just a matrix for data
        data.found <- FALSE
        if (any(c('rz','norm','log','transformed','raw') %in% tolower(names(ds$data)))) {
            data.found <- TRUE
        }

        if (!data.found) {
            message("ERROR   : 'rz','norm','log','transformed', OR 'raw' NOT found in data, must be ONE of them")
        }

        data_list <- get('data', ds)
        data_names <- ls(data_list)

        for (i in 1:length(data_names)) {
            cat('checking ', data_names[i])
            data_to_check <- paste0('data$', data_names[i])
            temp_data <- get(data_names[i], data_list)
            if (!is.numeric(temp_data)) {
                message(paste0(data_to_check, ' is not numeric'))
            }

            if (num_annot_samples != nrow(temp_data)) {
                cat('WARNING : number of samples (', num_annot_samples, ') != ',
                    'number of data["', data_names[i], '"] rows (', nrow(temp_data), ')\n')
            }

            if (NCOL(temp_data) != NROW(annot_ids)) {
                cat('WARNING : number of annotations (', NROW(annot_ids), ') != ',
                               'number of data["', data_names[i], '"] cols (', NCOL(temp_data), ')\n')

                x <- setdiff(annot_ids, colnames(temp_data))
                y <- setdiff(colnames(temp_data), annot_ids)

                if (length(x) > 0) {
                    #cat("WARNING : annotations with no data:", paste(x, sep = ",", collapse = ","), "\n")
                    cat("WARNING : # annotations with no data:", length(x), "\n")
                }

                if (length(y) > 0) {
                    #cat("WARNING : data with no annotations:", paste(y, sep = ",", collapse = ","), "\n")
                    cat("WARNING : # data with no annotations:", length(y), "\n")
                }
            }


            #perc_missing = (sum(is.na(ds_orig$data[data_names[i]])) * 100) / prod(dim(ds_orig$data[data_names[i]]))
            #cat("STATUS  : Percentage of missing original data: ", perc_missing, "\n")

            perc_missing = (sum(is.na(temp_data)) * 100) / prod(dim(temp_data))
            cat("STATUS  : Percentage of missing synchronized data: ", perc_missing, "\n")
        }
    }
}
