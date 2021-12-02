
# #############################################################################
#
# Utility functions
#
# #############################################################################


`%not in%` <- function(x, table) match(x, table, nomatch = 0L) == 0L


#' Check value for validity and return it or a default
#'
#' @param value value to check
#' @param default default value to use if value is "invalid"
#'
#' @return `value` if it is valid, `default` otherwise
#'
#' @noRd
nvl <- function(value, default) {
    if (gtools::invalid(value)) {
        return(default)
    }

    value
}


#' Convert value to numeric if possible
#'
#' @param value value to check
#' @param default default value to use if value is "invalid"
#'
#' @return `value` if it is numeric, `default` otherwise
#'
#' @noRd
nvl_int <- function(value, default) {
    tryCatch(
        {
            n <- as.numeric(value)
            if ((n %% 1) == 0) {
                return(n)
            }
        },
        error = function(cond) {
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    default
}


# #############################################################################
#
# Data methods
#
# #############################################################################


#' Get the dataset by id (a string).  annot.samples can be used at the top level
#' to share sample annotations among other datasets.
#'
#' @param ds a string, either 'dataset.name' or just 'name'
#'
#' @return the dataset element
#'
#' @export
get_dataset <- function(ds) {
    dataset <- NULL

    if (!is.character(ds)) {
        stop(sprintf("ds should be a string, but it's a '%s'", class(ds)))
    }

    if (exists(ds)) {
        dataset <- get(ds)
    } else {
        expanded <- sprintf("dataset.%s", ds)
        if (exists(expanded)) {
            dataset <- get(expanded)
        }
    }

    if (gtools::invalid(dataset)) {
        stop(sprintf("'%s' does not exist", ds))
    }

    samples <- NULL

    if (exists('annot.samples', dataset)) {
        samples <- dataset$annot.samples
    } else if (exists('annot.samples')) {
        samples <- get('annot.samples')
    }

    if (gtools::invalid(samples)) {
        stop(stop("Unable to find annot.samples"))
    }

    dataset$annot.samples <- samples

    return(dataset)
}

#' Check dataset to see if the datatype value is "phenotype"
#'
#' @param ds a dataset object
#'
#' @return `TRUE` if the datatype is phenotype, `FALSE` otherwise
#' @export
is_phenotype <- function(ds) {
    if ("datatype" %in% names(ds)) {
        if (startsWith(tolower(ds$datatype), "pheno")) {
            return(TRUE)
        }
    } else {
        message("datatype not found in dataset")
    }

    FALSE
}



#' Get the data from a data from the dataset.
#'
#' The data element in the dataset list should either be a matrix or a named
#' list with each eleemnt being a matrix.
#'
#' If `data_name` is not specified, the order that will be returned is:
#'     data (if a matrix)
#'     rz, norm, raw, log, transformed (if a named list).
#'
#' @param ds a datase object
#' @param data_name A string containing which data element from the dataset's
#'     data element.
#'
#' @return the data element
#' @export
get_data <- function(ds, data_name = NULL) {
    ret <- NULL

    # check if data is a matrix, but something else was requested
    if (is.matrix(ds$data) && (!is.null(data_name))) {
        stop(sprintf("Specified data '%s' not found in dataset", data_name))
    }

    if (!gtools::invalid(data_name)) {
        if (is.matrix(ds$data)) {
            stop(sprintf("Specified data '%s' not found in dataset", data_name))
        }

        # return the requested data element in the named list
        ret <- ds$data[[data_name]]
    } else {
        # Order of return is data (if matrix), than:
        # rz, norm, raw, log, transformed
        if (is.matrix(ds$data)) {
            ret <- ds$data
        } else {
            if (!is.null(ds$data$rz)) {
                ret <- ds$data$rz
            } else if (!is.null(ds$data$norm)) {
                ret <- ds$data$norm
            } else if (!is.null(ds$data$raw)) {
                ret <- ds$data$raw
            } else if (!is.null(ds$data$log)) {
                ret <- ds$data$log
            } else if (!is.null(ds$data$transformed)) {
                ret <- ds$data$transformed
            }
        }
    }

    if (gtools::invalid(ret)) {
        stop("Unable to find data in dataset")
    }

    ret
}


#' Use the existing covar.matrix or create it.
#'
#' @param ds the dataset object
#' @param id the phenotype identifier, for phenotype datasets
#'
#' @return the covar matrix
#'
#' @importFrom rlang .data
get_covar_matrix <- function(ds, id = NULL) {
    if (exists('covar.matrix', ds)) {
        covar <- ds$covar.matrix
    } else {
        # we can generate covar.matrix if it doesn't exist
        if (is_phenotype(ds)) {
            # get the annot.pheno row to get use.covar variable from the
            # annotations
            pheno <- ds$annot.pheno %>% dplyr::filter(.data$data_name == id)

            if (gtools::invalid(pheno)) {
                stop(sprintf("Cannot find phenotype '%s' in dataset", id))
            }

            # create a string (model formula) from the use.covar column
            formula_str <- paste0("~", gsub(":", "+", pheno$use_covar))
        } else {
            formula_str <- paste0(ds$covar.info$sample_column, collapse="+")
            formula_str <- paste0("~", formula_str)
        }

        # get the sample id field
        sample_id_field <- get_sample_id_field(ds)

        # convert samples to data.frame because QTL2 relies heavily
        # on rownames and colnames, rownames currently are or will
        # soon be deprecated in tibbles
        samples <- as.data.frame(ds$annot.samples)

        # set the rownames so scan1 will work
        rownames(samples) <-
            (samples %>%  dplyr::select(dplyr::matches(sample_id_field)))[[1]]

        # [, -1, drop = FALSE] will drop the (Intercept) column
        covar <- stats::model.matrix.lm(
            stats::as.formula(formula_str),
            data = samples,
            na.action = stats::na.pass
        )

        covar <- covar[, -1, drop = FALSE]
    }

    covar
}

#' Get all "dataset.*" information
#'
#' This will return a named list of all datasets and the ensmebl version.
#'
#' @return A named list of all the dataset objects, along with the
#'   ensembl.version.
#' @export
get_dataset_info <- function() {
    datasets <- grep('^dataset*', utils::apropos('dataset\\.'), value = TRUE)
    ret <- c()

    ensembl_version <-
        utils::apropos("^ensembl(.|_)version$", ignore.case = TRUE)

    if (length(ensembl_version) != 0) {
        ensembl_version <- get(ensembl_version)
    } else {
        ensembl_version <- NULL
    }

    for (d in datasets) {
        ds <- get(d)

        annotations <- list()

        if (tolower(ds$datatype) == 'mrna') {
            annot_mrna <- ds$annot.mrna %>% janitor::clean_names()
            annotations <- list(ids = annot_mrna$gene_id)
        } else if(tolower(ds$datatype) == 'protein') {
            annot_protein <- ds$annot.protein %>% janitor::clean_names()
            annotations <-
                list(ids = tibble::tibble(protein_id = annot_protein$protein_id,
                                          gene_id    = annot_protein$gene_id))
        } else if(is_phenotype(ds)) {
            annotations <-
                ds$annot.phenotype %>%
                janitor::clean_names() %>%
                dplyr::filter(.data$omit == FALSE)
        }

        covar_info <- ds$covar.info %>% janitor::clean_names()

        ds_ensembl_version <- ensembl_version

        temp_ensembl <- grep(
            "^ensembl(.|_)version$",
            names(ds),
            ignore.case = TRUE,
            value = TRUE
        )

        if (!gtools::invalid(temp_ensembl)) {
            ds_ensembl_version <- ds[[temp_ensembl]]
        }

        temp <- list(
            id              = d,
            annotations     = annotations,
            display_name    = nvl(ds$display.name, d),
            datatype        = ds$datatype,
            covar_info      = covar_info,
            ensembl_version = ds_ensembl_version
        )

        ret <- c(ret, list(temp))
    }

    list(datasets        = ret,
         ensembl_version = nvl(ensembl_version, NULL))
}

#' Get all "dataset.*" statistics.
#'
#' This will return a named list of all datasets and some statistics.
#'
#' @return A named list of all the dataset objects.
#' @export
get_dataset_stats <- function() {
    datasets <- grep('^dataset*', utils::apropos('dataset\\.'), value = TRUE)
    ret <- c()

    for (d in datasets) {
        ds <- get(d)

        num_annotations <- NA

        if (tolower(ds$datatype) == 'mrna') {
            num_annotations <- NROW(ds$annot.mrna)
        } else if(tolower(ds$datatype) == 'protein') {
            num_annotations <- NROW(ds$annot.protein)
        } else if(is_phenotype(ds)) {
            num_annotations <- NROW(ds$annot.phenotype)
        }

        temp <- list(id              = d,
                     display_name    = nvl(ds$display.name, d),
                     datatype        = ds$datatype,
                     num_annotations = num_annotations,
                     num_samples     = NROW(ds$annot.samples))

        ret <- c(ret, list(temp))
    }

    ret
}


