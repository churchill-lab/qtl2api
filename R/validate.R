#' Validate the dataset to make sure qtl2api can use it.
#'
#' @param dataset_id the dataset_id as a string identifier
#' @param extensive TRUE to perform scans on random elements in the data
#'
#' @export
validate_dataset <- function(dataset_id, extensive = FALSE) {
    cat("\nChecking '", dataset_id, "' ...\n", sep = '')

    # get the dataset
    ds_orig <- get_dataset_by_id(dataset_id)
    ds <- synchronize_dataset(ds_orig)

    if (gtools::invalid(ds)) {
        message(paste0("ERROR: dataset not found: ", dataset_id))
        return()
    }

    # dataset should be a list
    if (!is.list(ds_orig)) {
        message(paste0("ERROR: dataset should be a list, but found: ", class(ds_orig)))
    }

    if ('datatype' %not in% names(ds_orig)) {
        message("ERROR: dataset must contain 'datatype'")
    }

    if ('display.name' %not in% names(ds_orig)) {
        message("ERROR: dataset must contain 'display.name'")
    }

    datatype = ds[['datatype']]
    is_mrna <- FALSE
    is_protein <- FALSE
    is_pheno <- FALSE

    if (tolower(datatype) == 'mrna') {
        is_mrna <- TRUE
    } else if (tolower(datatype) == 'protein') {
        is_protein <- TRUE
    } else if (is_phenotype(ds)) {
        is_pheno <- TRUE
    } else {
        message(paste0("ERROR: datatype is invalid: ", datatype))
    }

    validate_annotations(dataset_id)

    validate_samples(dataset_id)

    validate_covar_info(dataset_id)

    validate_lod_peaks(dataset_id)

    validate_data(dataset_id)

    if (extensive) {
        validate_dataset_extensive(dataset_id)
    }
}

#' Perform scans on random data in the dataset designated by dataset_id.
#'
#' @param dataset_id the dataset as a string identifier
#'
#' @export
validate_dataset_extensive <- function(dataset_id) {
    # get the dataset
    dataset <- NA

    if (exists(dataset_id)) {
        dataset <- get(dataset_id)
    }

    ds <- synchronize_dataset(dataset)


    datatype = ds[['datatype']]
    id <- get_random_id(ds)

    covar_temp <-
        ds$covar_info
        dplyr::filter(.data$primary == TRUE)

    intcovar <- covar_temp$sample_column

    chrom <- sample(names(K), 1)

    markers_temp <- markers %>%
        dplyr::sample_n(1)

    marker <- markers_temp$marker.id

    cat("Checking get_peaks_all\n")
    tryCatch(
        {
            temp <- get_lod_peaks_all(ds)
        },
        error = function(cond) {
            message(cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    cat(paste0("Checking get_lod_scan ", id, "\n"))
    tryCatch(
        {
            temp <- get_lod_scan(ds, id)
        },
        error = function(cond) {
            message(cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )


    cat(paste0("Checking get_lod_scan ", id, " intcovar: ", intcovar, "\n"))
    tryCatch(
        {
            temp <- get_lod_scan(ds, id, intcovar)
        },
        error = function(cond) {
            message(cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    cat(paste0("Checking get_expression ", id, "\n"))
    tryCatch(
        {
            temp <- get_expression(ds, id)
        },
        error = function(cond) {
            message(cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    cat(paste0("Checking get_founder_coefficients ", id, " chrom: ", chrom, "\n"))
    tryCatch(
        {
            temp <- get_founder_coefficients(ds, id, chrom)
        },
        error = function(cond) {
            message(cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    cat(paste0("Checking get_founder_coefficients ", id, " chrom: ", chrom, ", intcovar: ", intcovar, "\n"))
    tryCatch(
        {
            temp <- get_founder_coefficients(ds, id, chrom, intcovar = intcovar)
        },
        error = function(cond) {
            message(cond$message)
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    if (!is_phenotype(ds)) {
        cat(paste0("Checking get_mediation ", id, " marker: ", marker, "\n"))
        tryCatch(
            {
                temp <- get_mediation(ds, id, marker)
            },
            error = function(cond) {
                message(cond$message)
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
        message("No datasets found!")
    }

    if (!exists('ensembl.version')) {
        message("ensembl.version does not exist")
    }

    # Check genoprobs and K.
    if(length(genoprobs) != length(K)) {
        message(paste0("genoprobs (", length(genoprobs), ") and K (", length(K), ") do not have the same length."))
    } else {
        if(any(names(genoprobs) != names(K))) {
            message("names of genoprobs and K do not match.")
        }

        rownames.eq <- mapply(function(x, y) { all(rownames(x) == rownames(y)) }, genoprobs, K)
        if(any(rownames.eq == FALSE)) {
            message("sample IDs do not match between genoprobs and K.")
        }
    }

    # Check Marker IDs for genoprobs and map
    if(length(genoprobs) != length(map)) {
        message(paste0("genoprobs (", length(genoprobs), ") and map (", length(map), ") do not have the same length."))
    } else {
        rownames.eq <- mapply(function(x, y) { all(dimnames(x)[[3]] == names(y)) }, genoprobs, map)
        if(any(rownames.eq == FALSE)) {
            message("marker names do not match between genoprobs and map.")
        }
    }

    # Check dimensions of markers and map.
    map.length = sum(sapply(map, length))
    if(map.length != nrow(markers)) {
        message(paste("Number of rows in markers (", nrow(markers), ") does not equal the number of markers in map (", map.length, ")"))
    }

    for (ds in datasets) {
        validate_dataset(ds, extensive)
    }
}


#' Validate the annotations in the dataset to make sure qtl2api can use it.
#'
#' @param dataset_id the dataset_id as a string identifier
#'
#' @export
validate_annotations <- function(dataset_id) {
    cat("Checking annotations\n")

    # get the dataset
    ds_orig <- get_dataset_by_id(dataset_id)
    ds <- synchronize_dataset(ds_orig)

    annots <- NULL
    annots_orig <- NULL

    if (tolower(ds$datatype) == "mrna") {
        if ('annot.mrna' %not in% names(ds_orig)) {
            message("ERROR: annot.mrna not found in dataset")
        }

        if (!tibble::is_tibble(ds_orig$annot.mrna)) {
            message("ERROR: annot.mrna should be a tibble, but found: ", class(ds$annot.mrna))
        }

        annots <- ds$annot_mrna
        annots_orig <- ds_orig$annot.mrna

        if (any(duplicated(annots$gene_id))) {
            message("ERROR: There are duplicated gene identifiers in annot.mrna")
        }
    } else if (tolower(ds$datatype) == "protein") {
        if ('annot.protein' %not in% names(ds_orig)) {
            message("ERROR: annot.protein not found in dataset")
        }

        if (!tibble::is_tibble(ds_orig$annot.protein)) {
            message(paste0("ERROR: annot.protein should be a tibble, but found: ", class(ds$annot.protein)))
        }

        annots <- ds$annot_protein
        annots_orig <- ds_orig$annot.protein

        if (any(duplicated(annots$protein_id))) {
            message("ERROR: There are duplicated protein identifiers in annot.protein")
        }
    } else if (is_phenotype(ds)) {
        if ('annot.phenotype' %not in% names(ds_orig)) {
            message("ERROR: annot.phenotype not found in dataset")
        }

        if (!tibble::is_tibble(ds_orig$annot.pheno)) {
            message(paste0("ERROR: annot.phenotype should be a tibble, but found ", class(ds_orig$annot.phenotype)))
        }

        annots <- ds$annot_phenotype
        annots_orig <- ds_orig$annot.phenotype

        if (any(duplicated(annots$data_name))) {
            message("ERROR: There are duplicated data_name annotations in annot.phenotype")
        }
    }


    if (is_phenotype(ds)) {
        expected_names <- c('data_name', 'short_name', 'description', 'is_id',
                            'category', 'is_numeric', 'is_date', 'is_factor',
                            'factor_levels', 'is_pheno', 'omit', 'use_covar')

        for (n in expected_names) {
            if (n %not in% names(annots)) {
                message(paste0("ERROR:", n, "not found in annot_phenotype"))
            }
        }

        expected_logical <- c('is_id', 'is_pheno', 'is_numeric', 'is_date',
                              'is_factor', 'omit')

        for (n in expected_logical) {
            if(!is.logical(annots[[n]])) {
                message(paste0("ERROR: annot_phenotype$", n, " should be logical, not ", class(annots[[n]])))
            }
        }

        if(is.logical(annots$is_id)) {
            # one and only 1 ID
            theID <- annots[which(annots$is_id == TRUE),]$data_name

            if(length(theID) != 1) {
                message('ERROR: annot_phenotype$is_id should have 1 and only 1 row set to TRUE')
            }
        }
    } else {
        annot_name <- NULL

        if (tolower(ds$datatype) == "protein") {
            annot_name <- 'annot_protein'

            if ('protein_id' %not in% names(annots)) {
                message("ERROR: protein_id not found in annot_protein")
            }
        } else {
            annot_name <- 'annot_mrna'
        }

        expected_names <-  c('gene_id', 'symbol', 'chr', 'start', 'end',
                             'strand', 'middle', 'nearest_marker_id')
        for (n in expected_names) {
            if (n %not in% names(annots)) {
                message(paste0("ERROR:", n, "not found in", annot_name))
            }
        }

        if (any(annots$start > 10000.0)) {
            message(paste0(annot_name, '$start should be in Mbp not bp', sep=''))
        } else if (any(annots$end > 10000.0)) {
            message(paste0(annot_name, '$end should be in Mbp not bp', sep=''))
        } else if (any(annots$middle > 10000.0)) {
            message(paste0(annot_name, '$middle should be in Mbp not bp', sep=''))
        }
    }

    num_annots_orig <- NROW(annots_orig)
    num_annots_synch <- NROW(annots)

    if (num_annots_orig != num_annots_synch) {
        cat("WARNING: synchronizing annotations, # samples changed from", num_annots_orig, "to", num_annots_synch, "\n")
    }
}


#' Validate the annot.samples in the dataset to make sure qtl2api can use it.
#'
#' @param dataset_id the dataset_id as a string identifier
#'
#' @export
validate_samples <- function(dataset_id) {
    cat("Checking samples\n")

    # get the dataset
    ds_orig <- get_dataset_by_id(dataset_id)
    ds <- synchronize_dataset(ds_orig)

    if ('annot.samples' %not in% names(ds_orig)) {
        message('ERROR: annot.samples not found')
    }

    if (!tibble::is_tibble(ds_orig$annot.sample)) {
        message(paste0("ERROR: annot.samples should be a tibble, but found ", class(ds$annot.sample)))
    }

    sample_id_field <- ds$sample_id_field
    if (gtools::invalid(sample_id_field)) {
        message('ERROR: unable to determine sample id field')
    }

    if (any(duplicated(ds$annot.sample[sample_id_field]))) {
        message("ERROR: there are duplicated annotations in annot.samples")
    }

    num_annots_orig <- NROW(ds_orig$annot.samples)
    num_annots_synch <- NROW(ds$annot_samples)

    if (num_annots_orig != num_annots_synch) {
        cat("WARNING: synchronizing samples, # samples changed from", num_annots_orig, "to", num_annots_synch, "\n")
    }
}


#' Validate the covar.info in the dataset to make sure qtl2api can use it.
#'
#' @param dataset_id the dataset_id as a string identifier
#'
#' @export
validate_covar_info <- function(dataset_id) {
    cat("Checking covar.info\n")

    # get the dataset
    ds_orig <- get_dataset_by_id(dataset_id)
    ds <- synchronize_dataset(ds_orig)

    if ('covar.info' %not in% names(ds_orig)) {
        message('covar.info not found')
    }

    if (!tibble::is_tibble(ds_orig$covar.info)) {
        message(paste0("covar.info should be a tibble, but found '", class(ds_orig$covar.info), "'"))
    }

    if ('sample.column' %not in% names(ds_orig$covar.info)) {
        message('sample.column not found in covar.info')
    } else if ('display.name' %not in% names(ds_orig$covar.info)) {
        message('display.name not found in covar.info')
    } else if ('interactive' %not in% names(ds_orig$covar.info)) {
        message('interactive not found in covar.info')
    } else if ('primary' %not in% names(ds_orig$covar.info)) {
        message('primary not found in covar.info')
    } else if ('lod.peaks' %not in% names(ds_orig$covar.info)) {
        message('lod.peaks not found in covar.info')
    }

    if(!is.logical(ds_orig$covar.info$primary)) {
        message('covar.info$primary should be logical, not ', class(ds_orig$covar.info$primary))
    }

    if(!is.logical(ds_orig$covar.info$interactive)) {
        message('covar.info$interactive should be logical, not ', class(ds_orig$covar.info$interactive))
    }

    for(i in 1:nrow(ds$covar_info)) {
        row <- ds$covar_info[i, ]

        if (row$sample_column %not in% colnames(ds$annot_samples)) {
            message(paste0("covar.info$sample.column ('", row$sample_column, "') is not a column name in annot.samples"))
        }

        if (gtools::invalid(row$display_name)) {
            message("covar.info$display_name needs to have a value")
        }

        if (class(row$interactive) == 'logical') {
            if (row$interactive) {
                if (gtools::invalid(row$lod_peaks)) {
                    message("covar.info$interactive is TRUE, but covar.info$lod_peaks is NA")
                } else {
                    # check for existence of lod_peaks
                    if (gtools::invalid(ds_orig$lod.peaks[[row$lod_peaks]])) {
                        message(paste0("covar.info$interactive is TRUE, but covar.info$lod_peaks ('", row$lod_peaks, "') is not in lod.peaks"))
                    }
                }
            } else {
                if (!gtools::invalid(row$lod_peaks)) {
                    message(paste0("covar.info$interactive is FALSE, but covar.info$lod_peaks ('", row$lod_peaks, "') is set"))
                }
            }
        } else {
            message("covar.info$interactive is INVALID, please set correctly")
        }
    }

    if (class(ds$covar_info$primary) == "logical") {
        if (!any(ds$covar_info$primary)) {
            message("covar.info$primary needs 1 value set to TRUE")
        }
    } else {
        message("covar.info$primary class is not logical (TRUE or FALSE), please set correctly")
    }

    cat("Checking covar.matrix\n")

    if (is_phenotype(ds)) {
        phenos <- ds$annot_phenotype %>% dplyr::filter(.data$is_pheno == TRUE)
        phenos <- phenos$data_name

        i <- 0
        for (x in phenos) {
            i <- i + 1
            tryCatch(
                {
                    if (i %% 500 == 0) {
                        cat("... completed", i, "checks out of", length(phenos), "\n")
                    }
                    temp <- get_covar_matrix(ds, x)
                },
                error = function(cond) {
                    message('Check phenotype ', x)
                    message('Unable to generate covar.matrix, check covar.info')
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
                message('Unable to generate covar.matrix, check covar.info')
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
#' @param dataset_id the dataset_id as a string identifier
#'
#' @export
validate_lod_peaks <- function(dataset_id) {
    cat("Checking lod.peaks\n")

    # get the dataset
    ds_orig <- get_dataset_by_id(dataset_id)
    ds <- synchronize_dataset(ds_orig)

    if ('lod.peaks' %not in% names(ds_orig)) {
        message('lod.peaks not found')
    }

    if (!is.list(ds_orig$lod.peaks)) {
        message(paste0("lod.peaks should be a list, but found ", class(ds_orig$lod.peaks)))
    }

    if ('additive' %not in% names(ds_orig$lod.peaks)) {
        message("additive should be an element in lod.peaks")
    }

    # get the rest
    for (i in seq(nrow(ds$covar_info))) {
        cov_inf <- ds$covar_info[i, ]

        # only look at interactive peaks
        if (cov_inf$interactive) {
            peaks <- ds_orig$lod.peaks[[cov_inf$lod_peaks]]

            if (gtools::invalid(peaks)) {
                message(paste0("Unable to find lod.peaks '", cov_inf$lod_peaks, "'"))
            } else {

                peaks <- peaks %>%
                    janitor::clean_names()


                if (tolower(ds$datatype) == "protein") {
                    if (length(setdiff(peaks$protein_id, ds$annot_protein$protein_id))) {
                        message(paste0('not all lod.peaks$protein_id are in annot.protein$protein_id'))
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        message(paste0('not all lod.peaks$marker_id are in markers'))
                    }
                } else if (tolower(ds$datatype) == "mrna") {
                    if (length(setdiff(peaks$gene_id, ds$annot_mrna$gene_id))) {
                        message(paste0('not all lod.peaks$gene_id are in annot.mrna$gene_id'))
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        message(paste0('not all lod.peaks$marker_id are in markers'))
                    }
                } else {
                    if (length(setdiff(peaks$data_name, ds$annot_phenotype$data_name))) {
                        message(paste0('not all lod.peaks["', cov_inf$sample_column, '"]$data_name are in annot.phenotype$data_name'))
                    }
                }
            }
        }
    }
}

#' Validate the data in the dataset to make sure qtl2api can use it.
#'
#' @param dataset_id the dataset_id as a string identifier
#'
#' @export
validate_data <- function(dataset_id) {
    cat("Checking data\n")

    # get the dataset
    ds_orig <- get_dataset_by_id(dataset_id)
    ds <- synchronize_dataset(ds_orig)

    num_annot_samples <- nrow(ds_orig$annot.samples)

    if (tolower(ds$datatype) == "mrna") {
        annot_ids <- ds$annot_mrna$gene_id
    } else if (tolower(ds$datatype) == "protein") {
        annot_ids <- ds$annot_protein$protein_id
    } else {
        annots_temp <-
            ds$annot_phenotype %>%
            dplyr::filter(
                .data$omit == FALSE,
                .data$is_pheno == TRUE
            )
        annot_ids <- annots_temp$data_name
    }

    if ('data' %not in% names(ds)) {
        message('data not found')
    }

    if (is.matrix(ds$data)) {
        # check if the data is numeric
        if (!is.numeric(ds$data)) {
            message(paste0(dataset_id, '$data matrix is not numeric'))
        }

        if (num_annot_samples != nrow(ds$data)) {
            cat(paste0('WRANING: number of samples (', num_annot_samples, ') != ',
                           'number of data rows (', nrow(ds$data), ')'))
        }


        if (NCOL(ds$data) != NROW(annot_ids)) {
            cat(paste0('WARNING: number of annotations (', NROW(annot_ids), ') != ',
                           'number of data cols (', NCOL(ds$data), ')'))

            x <- setdiff(annot_ids, colnames(ds$data))
            y <- setdiff(colnames(ds$data), annot_ids)

            if (length(x) > 0) {
                cat("WARNING: annotations with no data:")
                cat(paste(x, sep = ",", collapse = ","))
            }

            if (length(y) > 0) {
                cat("WARNING: data with no annotations:")
                cat(paste(y, sep = ",", collapse = ","))
            }
        }


    } else if (is.list(ds$data)) {
        data.found <- FALSE
        if (any(c('rz','norm','log','transformed','raw') %in% tolower(names(ds$data)))) {
            data.found <- TRUE
        }

        if (!data.found) {
            message("ERROR: 'rz','norm','log','transformed', OR 'raw' NOT found in data, must be ONE of them")
        }

        data_list <- get('data', ds)
        data_names <- ls(data_list)

        for (i in 1:length(data_names)) {
            data_to_check <- paste0(dataset_id, '$data$', data_names[i])
            temp_data <- get(data_names[i], data_list)
            if (!is.numeric(temp_data)) {
                message(paste0(data_to_check, ' is not numeric'))
            }

            if (num_annot_samples != nrow(temp_data)) {
                cat(paste0('number of samples (', num_annot_samples, ') != ',
                               'number of data["', data_names[i], '"] rows (', nrow(temp_data), ')'))
            }

            if (NCOL(temp_data) != NROW(annot_ids)) {
                cat(paste0('number of annotations (', NROW(annot_ids), ') != ',
                               'number of data["', data_names[i], '"] cols (', NCOL(temp_data), ')'))

                x <- setdiff(annot_ids, colnames(temp_data))
                y <- setdiff(colnames(temp_data), annot_ids)

                if (length(x) > 0) {
                    cat("WARNING: annotations with no data:")
                    cat(paste(x, sep = ",", collapse = ","))
                }

                if (length(y) > 0) {
                    cat("data with no annotations:")
                    cat(paste(y, sep = ",", collapse = ","))
                }

            }


            rm(temp_data)
        }
    }
}
