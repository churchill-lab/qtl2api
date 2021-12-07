#' Validate the dataset to make sure qtl2api can use it.
#'
#' @param dataset the dataset as a string identifier
#' @param extensive TRUE to perform scans on random elements in the data
#'
#' @export
validate_dataset <- function(dataset, extensive = FALSE) {
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

    cat("\nChecking '", dataset, "' ...\n", sep = '')

    # get the dataset
    ds_orig <- get_dataset_by_id(dataset)
    ds <- synchronize_dataset(ds_orig)

    if (gtools::invalid(ds)) {
        stop(paste0("dataset not found '", dataset, '"'))
    }

    if (!is.list(ds_orig)) {
        message(paste0("dataset should be a list, but found '", class(ds_orig), "'"))
    }

    if ('datatype' %not in% names(ds_orig)) {
        message("dataset must contain 'datatype'")
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
        message(paste0("datatype is '", datatype, "', but should be mRNA, protein, or phenotype"))
    }

    ###########################################################################
    #
    # annotations
    #
    ###########################################################################
    cat("Checking annotations\n")

    annots <- NULL

    if (is_mrna) {
        if ('annot.mrna' %not in% names(ds_orig)) {
            message(paste0("annot.mrna not found in '", dataset, "'"))
        }

        if (!tibble::is_tibble(ds_orig$annot.mrna)) {
            message(paste0("annot.mrna should be a tibble, but found '", class(ds$annot.mrna), "'"))
        }

        annots <- ds$annot_mrna

        if (any(duplicated(annots$gene_id))) {
            message("There are duplicated gene_id annotations in annot.mrna")
        }
    } else if (is_protein) {
        if ('annot.protein' %not in% names(ds_orig)) {
            message(paste0("annot.protein not found in '", dataset, "'"))
        }

        if (!tibble::is_tibble(ds_orig$annot.protein)) {
            message(paste0("annot.protein should be a tibble, but found ", class(ds$annot.protein)))
        }

        annots <- ds$annot_protein

        if (any(duplicated(annots$protein_id))) {
            message("There are duplicated protein_id annotations in annot.protein")
        }
    } else if (is_pheno) {
        if ('annot.phenotype' %not in% names(ds_orig)) {
            message(paste0("annot.phenotype not found in '", dataset, "'"))
        }

        if (!tibble::is_tibble(ds_orig$annot.pheno)) {
            message(paste0("annot.phenotype should be a tibble, but found ", class(ds$annot.phenotype)))
        }

        annots <- ds$annot_phenotype

        if (any(duplicated(annots$data_name))) {
            message("There are duplicated data_name annotations in annot.phenotype")
        }
    }

    if (is_pheno) {
        if ('data_name' %not in% names(annots)) {
            message('data_name not found in annot_phenotype')
        }
        if ('short_name' %not in% names(annots)) {
            message('short_name not found in annot_phenotype')
        }
        if ('description' %not in% names(annots)) {
            message('description not found in annot_phenotype')
        }
        if ('is_id' %not in% names(annots)) {
            message('is_id not found in annot_phenotype')
        }
        if ('category' %not in% names(annots)) {
            message('category not found in annot_phenotype')
        }
        if ('is_numeric' %not in% names(annots)) {
            message('is_numeric not found in annot_phenotype')
        }
        if ('is_date' %not in% names(annots)) {
            message('is_date not found in annot_phenotype')
        }
        if ('is_factor' %not in% names(annots)) {
            message('is_factor not found in annot_phenotype')
        }
        if ('factor_levels' %not in% names(annots)) {
            message('factor_levels not found in annot_phenotype')
        }
        if ('is_pheno' %not in% names(annots)) {
            message('is_pheno not found in annot_phenotype')
        }
        if ('omit' %not in% names(annots)) {
            message('omit not found in annot_phenotype')
        }
        if ('use_covar' %not in% names(annots)) {
            message('use_covar not found in annot_phenotype')
        }

        if(!is.logical(annots$is_id)) {
            message('annot_phenotype$is_id should be logical, not ', class(annots$is_id))
        }
        if(!is.logical(annots$is_pheno)) {
            message('annot_phenotype$is_pheno should be logical, not ', class(annots$is_pheno))
        }
        if(!is.logical(annots$is_numeric)) {
            message('annot_phenotype$is_numeric should be logical, not ', class(annots$is_numeric))
        }
        if(!is.logical(annots$is_date)) {
            message('annot_phenotype$is_date should be logical, not ', class(annots$is_date))
        }
        if(!is.logical(annots$is_factor)) {
            message('annot_phenotype$is_factor should be logical, not ', class(annots$is_factor))
        }
        if(!is.logical(annots$omit)) {
            message('annot_phenotype$omit should be logical, not ', class(annots$omit))
        }

        # one and only 1 ID
        theID <- annots[which(annots$is_id == TRUE),]$data_name

        if(length(theID) != 1) {
            message('annot_phenotype$is_id should have 1 and only 1 row set to TRUE')
        }

    } else {
        annot_name <- NULL

        if (is_protein) {
            annot_name <- 'annot_protein'

            if ('protein_id' %not in% names(annots)) {
                message('protein_id not found in annot_protein')
            }
        } else {
            annot_name <- 'annot_mrna'
        }

        if ('gene_id' %not in% names(annots)) {
            message(paste('gene_id not found in ', annot_name))
        } else if ('symbol' %not in% names(annots)) {
            message(paste('symbol not found in ', annot_name))
        } else if ('chr' %not in% names(annots)) {
            message(paste('chr not found in ', annot_name))
        } else if ('start' %not in% names(annots)) {
            message(paste('start not found in ', annot_name))
        } else if ('end' %not in% names(annots)) {
            message(paste('end not found in ', annot_name))
        } else if ('strand' %not in% names(annots)) {
            message(paste('strand not found in ', annot_name))
        } else if ('middle' %not in% names(annots)) {
            message(paste('strand not found in ', annot_name))
        } else if ('nearest_marker_id' %not in% names(annots)) {
            message(paste('nearest_marker_id not found in ', annot_name))
        }

        if (any(annots$start > 10000.0)) {
            message(paste0(annot_name, '$start should be in Mbp not bp', sep=''))
        } else if (any(annots$end > 10000.0)) {
            message(paste0(annot_name, '$end should be in Mbp not bp', sep=''))
        } else if (any(annots$middle > 10000.0)) {
            message(paste0(annot_name, '$middle should be in Mbp not bp', sep=''))
        }
    }

    ###########################################################################
    #
    # samples
    #
    ###########################################################################
    cat("Checking annot.samples\n")

    print(names(ds))
    if ('annot.samples' %not in% names(ds_orig)) {
        message('annot.samples not found')
    }

    if (!tibble::is_tibble(ds_orig$annot.sample)) {
        message(paste0("annot.samples should be a tibble, but found ", class(ds$annot.sample)))
    }

    sample_id_field <- get_sample_id_field(ds_orig)
    if (gtools::invalid(sample_id_field)) {
        message('unable to determine sample id field')
    }

    if (any(duplicated(ds$annot.sample[sample_id_field]))) {
        message("There are duplicated annotations in annot.samples")
    }

    ###########################################################################
    #
    # display.name
    #
    ###########################################################################
    cat("Checking display.name\n")

    if ('display.name' %not in% names(ds_orig)) {
        message('display.name not found')
    }

    ###########################################################################
    #
    # covar.info
    #
    ###########################################################################
    cat("Checking covar.info\n")

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

    for(i in 1:nrow(covar_info)) {
        row <- covar_info[i, ]

        if (row$sample_column %not in% colnames(ds$annot.samples)) {
            message(paste0("covar.info$sample_column ('", row$sample_column, "') is not a column name in annot.samples"))
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
                    if (gtools::invalid(ds$lod.peaks[[row$lod_peaks]])) {
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

    if (class(covar_info$primary) == "logical") {
        if (!any(covar_info$primary)) {
            message("covar.info$primary needs 1 value set to TRUE")
        }
    } else {
        message("covar.info$primary class is not logical (TRUE or FALSE), please set correctly")
    }


    ###########################################################################
    #
    # covar.matrix
    #
    ###########################################################################
    cat("Checking covar.matrix\n")

    if (is_phenotype(ds)) {
        phenos <- annots %>% dplyr::filter(.data$is_pheno == TRUE)
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


    ###########################################################################
    #
    # lod.peaks
    #
    ###########################################################################
    cat("Checking lod.peaks\n")

    if ('lod.peaks' %not in% names(ds)) {
        message('lod.peaks not found')
    }

    if (!is.list(ds$lod.peaks)) {
        message(paste0("lod.peaks should be a list, but found ", class(ds$lod.peaks)))
    }

    if ('additive' %not in% names(ds$lod.peaks)) {
        message("additive should be an element in lod.peaks")
    }

    #
    # lod.peaks
    #

    # get the rest
    for (i in seq(nrow(covar_info))) {
        cov_inf <- covar_info[i, ]

        # only look at interactive peaks
        if (cov_inf$interactive) {
            peaks <- ds$lod.peaks[[cov_inf$lod_peaks]]

            if (gtools::invalid(peaks)) {
                message(paste0("Unable to find lod.peaks '", cov_inf$lod_peaks, "'"))
            } else {

                peaks <- peaks %>%
                    janitor::clean_names()


                if (is_pheno) {
                    if (length(setdiff(peaks$data_name, annots$data_name))) {
                        message(paste0('not all lod.peaks["', cov_inf$sample_column, '"]$data_name are in annot.phenotype$data_name'))
                    }
                } else if (is_protein) {
                    if (length(setdiff(peaks$protein_id, annots$protein_id))) {
                        message(paste0('not all lod.peaks$protein_id are in annot.protein$protein_id'))
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        message(paste0('not all lod.peaks$marker_id are in markers'))
                    }
                } else if (is_mrna) {
                    if (length(setdiff(peaks$gene_id, annots$gene_id))) {
                        message(paste0('not all lod.peaks$gene_id are in annot.mrna$gene_id'))
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        message(paste0('not all lod.peaks$marker_id are in markers'))
                    }
                }
            }
        }
    }

    ###########################################################################
    #
    # data
    #
    ###########################################################################
    cat("Checking data\n")

    num_annot_samples <- nrow(ds$annot.samples)

    if (is_mrna) {
        annot_ids <- annots$gene_id
    } else if (is_protein) {
        anno_ids <- annots$protein_id
    } else {
        annots_temp <-
            annots %>%
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
            message(paste0(dataset, '$data matrix is not numeric'))
        }

        if (num_annot_samples != nrow(ds$data)) {
            message(paste0('number of samples (', num_annot_samples, ') != ',
                           'number of data rows (', nrow(ds$data), ')'))
        }


        if (NCOL(ds$data) != NROW(annot_ids)) {
            message(paste0('number of annotations (', NROW(annot_ids), ') != ',
                           'number of data cols (', NCOL(ds$data), ')'))

            x <- setdiff(annot_ids, colnames(ds$data))
            y <- setdiff(colnames(ds$data), annot_ids)

            if (length(x) > 0) {
                message("annotations with no data:")
                message(paste(x, sep = ",", collapse = ","))
            }

            if (length(y) > 0) {
                message("data with no annotations:")
                message(paste(y, sep = ",", collapse = ","))
            }
        }


    } else if (is.list(ds$data)) {
        data.found <- FALSE
        if (any(c('rz','norm','log','transformed','raw') %in% tolower(names(ds$data)))) {
            data.found <- TRUE
        }

        if (!data.found) {
            message("'rz','norm','log','transformed', OR 'raw' NOT found in data, must be ONE of them")
        }

        data_list <- get('data', ds)
        data_names <- ls(data_list)

        for (i in 1:length(data_names)) {
            data_to_check <- paste0(dataset, '$data$', data_names[i])
            temp_data <- get(data_names[i], data_list)
            if (!is.numeric(temp_data)) {
                message(paste0(data_to_check, ' is not numeric'))
            }

            if (num_annot_samples != nrow(temp_data)) {
                message(paste0('number of samples (', num_annot_samples, ') != ',
                               'number of data["', data_names[i], '"] rows (', nrow(temp_data), ')'))
            }

            if (NCOL(temp_data) != NROW(annot_ids)) {
                message(paste0('number of annotations (', NROW(annot_ids), ') != ',
                               'number of data["', data_names[i], '"] cols (', NCOL(temp_data), ')'))

                x <- setdiff(annot_ids, colnames(temp_data))
                y <- setdiff(colnames(temp_data), annot_ids)

                if (length(x) > 0) {
                    message("annotations with no data:")
                    message(paste(x, sep = ",", collapse = ","))
                }

                if (length(y) > 0) {
                    message("data with no annotations:")
                    message(paste(y, sep = ",", collapse = ","))
                }

            }


            rm(temp_data)
        }
    }

    rm(num_annot_samples)

    if (extensive) {
        validate_dataset_extensive(dataset)
    }
}

#' Perform scans on random data in the datset.
#'
#' @param dataset the dataset as a string identifier
#'
#' @export
validate_dataset_extensive <- function(dataset) {
    # get the dataset
    ds <- NA

    if (exists(dataset)) {
        ds <- get(dataset)
    }

    datatype = ds[['datatype']]
    id <- get_random_id(ds)

    covar_temp <-
        ds$covar.info %>%
        janitor::clean_names() %>%
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


