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
#' @export
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
#' @export
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


#' Synchronize data and subset accordingly.
#'
#' @param dataset the dataset object
#'
#' @return list with 3 elements: `annots`, `samples` and `data`.
#'
#' @export
synchronize_data <- function(dataset) {
    # first thing is to get the annotation ids
    if (tolower(dataset$datatype) == 'mrna') {
        annots <-
            dataset$annot.mrna %>%
            janitor::clean_names()

        annot_ids <- annots$gene_id
    } else if (tolower(dataset$datatype) == 'protein') {
        annots <-
            dataset$annot.protein %>%
            janitor::clean_names()

        annot_ids <- annots$protein_id
    } else if (is_phenotype(dataset)) {
        annots <-
            dataset$annot.phenotype %>%
            janitor::clean_names() %>%
            dplyr::filter(.data$omit == FALSE, .data$is_pheno == TRUE)

        annot_ids <- annots$data_name
    } else {
        message(paste0("datatype is invalid: '", dataset$datatype, "'"))
    }

    # grab the data, row names are sample ids, column names are annotation ids
    data <- get_data(dataset)

    # get the sample id field and the intersecting sample ids
    sample_id_field <- get_sample_id_field(dataset)
    sample_ids <- intersect(
        rownames(data),
        dataset$annot.samples[[sample_id_field]]
    )

    if (length(sample_ids) == 0) {
        message("There are no samples in common")
    }

    # get the intersecting annotation ids
    annot_ids <- intersect(
        colnames(data),
        annot_ids
    )

    if (length(annot_ids) == 0) {
        message("There are no annotations in common")
    }

    # sort the ids to make sure the come back in order
    annot_ids <- sort(annot_ids)
    sample_ids <- sort(sample_ids)

    # filter the annots
    if (tolower(dataset$datatype) == 'mrna') {
        annots <- annots %>% dplyr::filter(.data$gene_id %in% annot_ids)
    } else if (tolower(dataset$datatype) == 'protein') {
        annots <- annots %>% dplyr::filter(.data$protein_id %in% annot_ids)
    } else if (is_phenotype(dataset)) {
        annots <- annots %>% dplyr::filter(.data$data_name %in% annot_ids)
    }

    # filter the data
    data <- data[sample_ids, annot_ids, drop = FALSE]

    # filter the samples
    samples <-
        dataset$annot.samples %>%
        dplyr::filter(!!as.name(sample_id_field) %in% sample_ids)

    list(
        annots  = annots,
        data    = data,
        samples = samples
    )
}


#' Get the dataset that is altered to be better use in the environment.
#'
#' @param dataset the dataset element
#'
#' @return the altered dataset element
#'
#' @export
synchronize_dataset <- function(dataset) {
    if (gtools::invalid(dataset)) {
        stop("invalid dataset object")
    } else if (class(dataset) != "list") {
        stop(paste0("dataset should be a list, but it's a ", class(dataset)))
    }

    if (!gtools::invalid(dataset$is_synchronized)) {
        # already synchronized
        return(dataset)
    }

    # TODO: samples can be top level, but not implemented

    # synchronize the data element
    ds_synch <- synchronize_data(dataset)

    # fix the covar.info names
    covar_info <-
        dataset$covar.info %>%
        janitor::clean_names()

    # fix the ensembl version
    ensembl_version <-
        utils::apropos("^ensembl(.|_)version$", ignore.case = TRUE)

    if (gtools::invalid(ensembl_version)) {
        ensembl_version <- get(ensembl_version)
    } else {
        ensembl_version <- NULL
    }

    ds_ensembl_version <- ensembl_version

    temp_ensembl <- grep(
        "^ensembl(.|_)version$",
        names(dataset),
        ignore.case = TRUE,
        value = TRUE
    )

    if (!gtools::invalid(temp_ensembl)) {
        ds_ensembl_version <- dataset[[temp_ensembl]]
    }

    sample_id_field = get_sample_id_field(dataset)

    ds <- list(
        annot_samples   = ds_synch$samples,
        covar_info      = covar_info,
        data            = ds_synch$data,
        datatype        = dataset$datatype,
        display_name    = dataset$display.name,
        ensembl_version = ensembl_version,
        sample_id_field = sample_id_field,
        is_synchronized = TRUE
    )

    if (tolower(dataset$datatype) == 'mrna') {
        ds$annot_mrna <- ds_synch$annots
    } else if (tolower(dataset$datatype) == 'protein') {
        ds$annot_protein <- ds_synch$annots
    } else if (startsWith(tolower(dataset$datatype), "pheno")) {
        ds$annot_phenotype <- ds_synch$annots

        ds$annot_phenotpe_extra <-
            dataset$annot.phenotype %>%
            janitor::clean_names() %>%
            dplyr::filter(.data$omit == FALSE & .data$is_pheno == FALSE)

    } else {
        message(paste0("datatype is invalid: '", dataset$datatype, "'"))
    }

    return(ds)
}


#' Get the dataset by id (a string).  annot.samples can be used at the top level
#' to share sample annotations among other datasets.
#'
#' @param dataset_id a string, either 'dataset.name' or just 'name'
#'
#' @return the dataset element
#'
#' @export
get_dataset_by_id <- function(dataset_id) {
    if (!is.character(dataset_id)) {
        stop(paste0("dataset_id should be a string, not ", class(dataset_id)))
    }

    dataset <- NULL

    if (exists(dataset_id)) {
        dataset <- get(dataset_id)
    } else {
        expanded <- sprintf("dataset.%s", dataset_id)
        if (exists(expanded)) {
            dataset <- get(expanded)
        }
    }

    if (gtools::invalid(dataset)) {
        stop(sprintf("'%s' does not exist", dataset_id))
    }

    return(dataset)
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
#' @param ds a dataset object
#' @param data_name A string containing which data element from the dataset's
#'     data element
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


#' Get a random identifier from the dataset.
#'
#' @param dataset a dataset object
#'
#' @return a random identifier
#' @export
get_random_id <- function(dataset) {
    if (tolower(dataset$datatype) == "mrna") {
        if ("annot.mrna" %in% dataset) {
            annot_ids <-
                dataset$annot.mrna %>%
                janitor::clean_names()
        } else {
            annot_ids <- dataset$annot_mrna
        }

        annot_ids <- annot_ids$gene_id
    } else if (tolower(dataset$datatype) == "protein") {
        if ("annot.protein" %in% dataset) {
            annot_ids <-
                dataset$annot.protein %>%
                janitor::clean_names()
        } else {
            annot_ids <- dataset$annot_protein
        }

        annot_ids <- annot_ids$protein_id
    } else if (is_phenotype(dataset)) {
        if ("annot.phenotype" %in% dataset) {
            annot_ids <-
                dataset$annot.phenotype %>%
                janitor::clean_names() %>%
                dplyr::filter(.data$is_pheno == TRUE)
        } else {
            annot_ids <- dataset$annot_phenotype
        }

        annot_ids <- annot_ids$data_name
    } else {
        stop("Unknown datatype: ", dataset$datatype)
    }

    data_ids <- colnames(get_data(dataset))

    # invalid_ids <- union(
    #     setdiff(data_ids, annot_ids),
    #     setdiff(annot_ids, data_ids)
    # )

    valid_ids <- intersect(annot_ids, data_ids)

    return(sample(valid_ids, 1))
}


#' Create the covar.matrix element.
#'
#' @param dataset the dataset object (synchronized)
#' @param id the phenotype identifier, for phenotype datasets
#'
#' @return the covar matrix
#'
#' @importFrom rlang .data
get_covar_matrix <- function(dataset, id = NULL) {
    if (is_phenotype(dataset)) {
        # get the annot.pheno row to get use.covar variable from the
        # annotations
        pheno <-
            dataset$annot_phenotype %>%
            dplyr::filter(.data$data_name == id)

        if (gtools::invalid(pheno)) {
            stop(sprintf("Cannot find phenotype '%s' in dataset", id))
        }

        # create a string (model formula) from the use.covar column
        formula_str <- paste0("~", gsub(":", "+", pheno$use_covar))
    } else {
        formula_str <- paste0(dataset$covar_info$sample_column, collapse="+")
        formula_str <- paste0("~", formula_str)
    }

    # get the sample id field
    sample_id_field <- dataset$sample_id_field

    # convert samples to data.frame because QTL2 relies heavily
    # on rownames and colnames, rownames currently are or will
    # soon be deprecated in tibbles
    samples <- as.data.frame(dataset$annot_samples)

    # set the rownames so scan1 will work
    rownames(samples) <-
        (samples %>%  dplyr::select(dplyr::matches(sample_id_field)))[[1]]

    # [, -1, drop = FALSE] will drop the (Intercept) column
    covar <- stats::model.matrix.lm(
        stats::as.formula(formula_str),
        data = samples,
        na.action = stats::na.pass
    )

    return(covar[, -1, drop = FALSE])
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

        # make sure samples and annotations are available
        ds_synchronized <- synchronize_data(ds)

        annotations <- list()

        if (tolower(ds$datatype) == 'mrna') {
            annotations <- list(ids = ds_synchronized$annots$gene_id)
        } else if(tolower(ds$datatype) == 'protein') {
            annotations <-
                list(
                    ids = tibble::tibble(
                        protein_id = ds_synchronized$annots$protein_id,
                        gene_id    = ds_synchronized$annots$gene_id
                    )
                )
        } else if(is_phenotype(ds)) {
            # this is trickier, we need to send back the is_pheno = FALSE too
            annotations <-
                ds$annot.phenotype %>%
                janitor::clean_names() %>%
                dplyr::filter(.data$omit == FALSE & .data$is_pheno == FALSE)

            annotations <- dplyr::bind_rows(
                annotations,
                ds_synchronized$annots
            )
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


