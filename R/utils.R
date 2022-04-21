# #############################################################################
#
# Utility functions
#
# #############################################################################

`%not in%` <- function(x, table) match(x, table, nomatch = 0L) == 0L


#' Get the version on qtl2api.
#'
#' @return version of qtl2api.
#'
#' @export
version <- function() {
    v <- unlist(utils::packageVersion("qtl2api"))
    paste(v, collapse=".")
}


#' Test if `val` is invalid.
#'
#' @param val value to be tested
#' @return TRUE if invalid, FALSE otherwise
#'
#' @export
invalid <- function(val) {
    return(gtools::invalid(val))
}


#' Test if `val` is valid.
#'
#' @param x value to be tested
#' @return TRUE if valid, FALSE otherwise
#'
#' @export
valid <- function(val) {
    return(!gtools::invalid(val))
}


#' Check value for validity and return it or a default
#'
#' @param value value to check
#' @param default default value to use if value is "invalid"
#'
#' @return `value` if it is valid, `default` otherwise
#'
#' @export
nvl <- function(value, default) {
    if (invalid(value)) {
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
        annots_field <- grep("^annots?(\\.|_){1}mrnas?$",
                             names(dataset),
                             value = TRUE)

        annots <-
            dataset[[annots_field]] %>%
            janitor::clean_names()

        annot_ids <- annots$gene_id
    } else if (tolower(dataset$datatype) == 'protein') {
        annots_field <- grep("^annots?(\\.|_){1}proteins?$",
                             names(dataset),
                             value = TRUE)

        annots <-
            dataset[[annots_field]] %>%
            janitor::clean_names()

        annot_ids <- annots$protein_id
    } else if (tolower(dataset$datatype) == 'phos') {
        annots_field <- grep("^annots?(\\.|_){1}phos?$",
                             names(dataset),
                             value = TRUE)

        annots <-
            dataset[[annots_field]] %>%
            janitor::clean_names()

        annot_ids <- annots$phos_id
    } else if (is_phenotype(dataset)) {
        annots_field <- grep("^annots?(\\.|_){1}pheno(type)?s?$",
                             names(dataset),
                             value = TRUE)

        annots <-
            dataset[[annots_field]] %>%
            janitor::clean_names() %>%
            dplyr::filter(.data$omit == FALSE, .data$is_pheno == TRUE)

        annot_ids <- annots$data_name
    } else {
        message(paste0("datatype is invalid: '", dataset$datatype, "'"))
    }

    # grab the data, row names are sample ids, column names are annotation ids
    data <- get_data(dataset)

    # get the sample id field and the intersecting sample ids
    annots_field <- grep("^annots?(\\.|_){1}samples?$",
                         names(dataset),
                         value = TRUE)

    sample_id_field <- get_sample_id_field(dataset)
    sample_ids <- intersect(
        rownames(data),
        dataset[[annots_field]][[sample_id_field]]
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
    } else if (tolower(dataset$datatype) == 'phos') {
        annots <- annots %>% dplyr::filter(.data$phos_id %in% annot_ids)
    } else if (is_phenotype(dataset)) {
        annots <- annots %>% dplyr::filter(.data$data_name %in% annot_ids)
    }

    # filter the data
    data <- data[sample_ids, annot_ids, drop = FALSE]

    # filter the samples
    annots_field <- grep("^annots?(\\.|_){1}samples?$",
                         names(dataset),
                         value = TRUE)

    samples <-
        dataset[[annots_field]] %>%
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
    if (invalid(dataset)) {
        stop("invalid dataset object")
    } else if (class(dataset) != "list") {
        stop(paste0("dataset should be a list, but it's a ", class(dataset)))
    }

    if (valid(dataset$is_synchronized)) {
        # already synchronized
        return(dataset)
    }

    # TODO: samples can be top level, but not implemented

    # synchronize the data element
    ds_synch <- synchronize_data(dataset)

    # fix the covar_info names
    covar_info <- NULL

    annots_field <- grep("^covar?(\\.|_){1}info$",
                         names(dataset),
                         value = TRUE)

    if ((length(annots_field) > 0) && (!is.null(dataset[[annots_field]]))) {
        covar_info <- dataset[[annots_field]]
    }

    if (!is.null(covar_info)) {
        covar_info <- covar_info %>% janitor::clean_names()
    }

    # fix the ensembl version
    ensembl_version <-
        utils::apropos("^ensembl(\\.|_){1}version$", ignore.case = TRUE)

    if (valid(ensembl_version)) {
        ensembl_version <- get(ensembl_version)
    } else {
        ensembl_version <- NULL
    }

    ds_ensembl_version <- ensembl_version

    temp_ensembl <- grep(
        "^ensembl(\\.|_){1}version$",
        names(dataset),
        ignore.case = TRUE,
        value = TRUE
    )

    if (valid(temp_ensembl)) {
        ds_ensembl_version <- dataset[[temp_ensembl]]
    }

    sample_id_field = get_sample_id_field(dataset)

    display_name_field <- grep(
        "^display(\\.|_){1}name$",
        names(dataset),
        ignore.case = TRUE,
        value = TRUE
    )

    datatype <- tolower(dataset$datatype)
    if (startsWith(datatype, "pheno")) {
        datatype <- 'phenotype'
    }

    ds <- list(
        annot_samples   = ds_synch$samples,
        covar_info      = covar_info,
        data            = ds_synch$data,
        datatype        = datatype,
        display_name    = dataset[[display_name_field]],
        ensembl_version = ensembl_version,
        sample_id_field = sample_id_field,
        is_synchronized = TRUE
    )

    if (tolower(dataset$datatype) == 'mrna') {
        if (all(ds_synch$annots$start < 1000)) {
            ds_synch$annots$start <- ds_synch$annots$start * 1000000
        }

        if (all(ds_synch$annots$end < 1000)) {
            ds_synch$annots$end <- ds_synch$annots$end * 1000000
        }

        ds$annot_mrna <- ds_synch$annots
    } else if (tolower(dataset$datatype) == 'protein') {
        if (all(ds_synch$annots$start < 1000)) {
            ds_synch$annots$start <- ds_synch$annots$start * 1000000
        }

        if (all(ds_synch$annots$end < 1000)) {
            ds_synch$annots$end <- ds_synch$annots$end * 1000000
        }

        ds$annot_protein <- ds_synch$annots
    } else if (tolower(dataset$datatype) == 'phos') {
        if (all(ds_synch$annots$start < 1000)) {
            ds_synch$annots$start <- ds_synch$annots$start * 1000000
        }

        if (all(ds_synch$annots$end < 1000)) {
            ds_synch$annots$end <- ds_synch$annots$end * 1000000
        }

        ds$annot_phos <- ds_synch$annots
    } else if (startsWith(tolower(dataset$datatype), "pheno")) {
        ds$annot_phenotype <- ds_synch$annots

        #ds$annot_phenotpe_extra <-
        #    dataset$annot_phenotype %>%
        #    janitor::clean_names() %>%
        #    dplyr::filter(.data$omit == FALSE & .data$is_pheno == FALSE)

    } else {
        message(paste0("datatype is invalid: '", dataset$datatype, "'"))
    }

    return(ds)
}


#' Get the dataset by id (a string).  annot_samples can be used at the top level
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

    if (invalid(dataset)) {
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
#' @param ds a dataset object (synchronized or not)
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

    if (valid(data_name)) {
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

    if (invalid(ret)) {
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
        if (invalid(dataset$is_synchronized)) {
            annots_field <- grep("^annots?(\\.|_){1}mrnas?$",
                                 names(dataset),
                                 value = TRUE)
            annot_ids <-
                dataset[[annots_field]] %>%
                janitor::clean_names()
        } else {
            annot_ids <- dataset$annot_mrna
        }

        annot_ids <- annot_ids$gene_id
    } else if (tolower(dataset$datatype) == "protein") {
        if (invalid(dataset$is_synchronized)) {
            annots_field <- grep("^annots?(\\.|_){1}proteins?$",
                                 names(dataset),
                                 value = TRUE)
            annot_ids <-
                dataset[[annots_field]] %>%
                janitor::clean_names()
        } else {
            annot_ids <- dataset$annot_protein
        }

        annot_ids <- annot_ids$protein_id
    } else if (tolower(dataset$datatype) == "phos") {
        if (invalid(dataset$is_synchronized)) {
            annots_field <- grep("^annots?(\\.|_){1}phos$",
                                 names(dataset),
                                 value = TRUE)
            annot_ids <-
                dataset[[annots_field]] %>%
                janitor::clean_names()
        } else {
            annot_ids <- dataset$annot_phos
        }

        annot_ids <- annot_ids$phos_id
    } else if (is_phenotype(dataset)) {
        if (invalid(dataset$is_synchronized)) {
            annots_field <- grep("^annots?(\\.|_){1}pheno(type)?s?$",
                                 names(dataset),
                                 value = TRUE)
            annot_ids <-
                dataset[[annots_field]] %>%
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
#' @return a named list list with 2 elements, covar_matrix and covar_formula
#'
#' @importFrom rlang .data
get_covar_matrix <- function(dataset, id = NULL) {
    # make sure samples and annotations are available
    ds <- synchronize_dataset(dataset)

    covar_formula <- NULL
    covar_matrix <- NULL

    if (is_phenotype(ds)) {
        # get the annot_phenotype row to get use_covar variable from the
        # annotations
        pheno <-
            ds$annot_phenotype %>%
            dplyr::filter(.data$data_name == id)

        if (invalid(pheno)) {
            stop(sprintf("Cannot find phenotype '%s' in dataset", id))
        }

        if (!is.null(pheno$use_covar)) {
            # create a string (model formula) from the use.covar column
            covar_formula <- paste0("~", gsub(":", "+", pheno$use_covar))
        }
    } else {
        if (!is.null(ds$covar_info)) {
            covar_formula <- paste0(ds$covar_info$sample_column, collapse="+")
            covar_formula <- paste0("~", covar_formula)
        }
    }

    if (!is.null(covar_formula)) {
        # get the sample id field
        sample_id_field <- ds$sample_id_field

        # convert samples to data.frame because QTL2 relies heavily
        # on rownames and colnames, rownames currently are or will
        # soon be deprecated in tibbles
        samples <- as.data.frame(ds$annot_samples)

        # create the model matrix, we use na.action = stats::na.pass so we can set
        # the rownames below.  We than use na.omit to filter down the data.
        covar_matrix <- stats::model.matrix.lm(
            stats::as.formula(covar_formula),
            data = samples,
            na.action = stats::na.pass
        )

        # drop the Intercept column
        covar_matrix <- covar_matrix[, -1, drop = FALSE]

        # drop the covar column if it has all identical values
        covar_matrix <- covar_matrix %>%
            tibble::as_tibble() %>%
            dplyr::select_if(function(col) length(unique(col))>1)

        # convert to a matrix and set the rownames so scan1 will work
        covar_matrix <- as.matrix(covar_matrix)

        rownames(covar_matrix) <-
            (samples %>% dplyr::select(dplyr::matches(sample_id_field)))[[1]]

        # do not need NA values
        covar_matrix <- stats::na.omit(covar_matrix)
    }

    list(
        covar_formula  = covar_formula,
        covar_matrix   = covar_matrix
    )
}


#' Get all "dataset.*" information
#'
#' This will return a named list of all datasets and the ensmebl version.
#'
#' @return A named list of all the dataset objects, along with the
#'   ensembl.version.
#' @export
get_dataset_info <- function() {
    datasets <- utils::apropos('^dataset\\.*', ignore.case = TRUE)
    ret <- c()

    ensembl_version_field <-
        utils::apropos("^ensembl(\\.|_){1}version$", ignore.case = TRUE)

    if (length(ensembl_version_field) != 0) {
        ensembl_version <- get(ensembl_version_field)
    } else {
        ensembl_version <- NULL
    }

    for (d in datasets) {
        ds <- get(d)

        # make sure samples and annotations are available
        ds_synchronized <- synchronize_data(ds)

        annotations <- list()

        if (ds$datatype == 'mrna') {
            annotations <-
                tibble::tibble(
                    gene_id    = ds_synchronized$annots$gene_id
                )
        } else if(ds$datatype == 'protein') {
            annotations <-
                tibble::tibble(
                    protein_id = ds_synchronized$annots$protein_id,
                    gene_id    = ds_synchronized$annots$gene_id
                )
        } else if(ds$datatype == 'phos') {
            annotations <-
                tibble::tibble(
                    phos_id    = ds_synchronized$annots$phos_id,
                    protein_id = ds_synchronized$annots$protein_id,
                    gene_id    = ds_synchronized$annots$gene_id
                )
        } else if(ds$datatype == 'phenotype') {
            # this is trickier, we need to send back the is_pheno = FALSE too
            # TODO: Rethink this, do we need is_pheno == FALSE?
            #
            # annotations <-
            #     ds$annot_phenotype %>%
            #     janitor::clean_names() %>%
            #     dplyr::filter(.data$omit == FALSE & .data$is_pheno == FALSE)
            #
            # annotations <- dplyr::bind_rows(
            #     annotations,
            #     ds_synchronized$annots
            # )
            annotations <- ds_synchronized$annots
        }

        annots_field <- grep("^covar?(\\.|_){1}info$",
                             names(ds),
                             value = TRUE)

        covar_info <- NULL

        if ((length(annots_field) != 0) && (!is.null(ds[[annots_field]]))) {
            covar_info <- ds[[annots_field]] %>% janitor::clean_names()
        }

        display_name_field <- grep(
            "^display(\\.|_){1}name$",
            names(ds),
            ignore.case = TRUE,
            value = TRUE
        )

        display_name <- d

        if (length(display_name_field) != 0) {
            display_name <- ds[[display_name_field]]
        }

        ds_ensembl_version <- ensembl_version

        temp_ensembl <- grep(
            "^ensembl(\\.|_){1}version$",
            names(ds),
            ignore.case = TRUE,
            value = TRUE
        )

        if (valid(temp_ensembl)) {
            ds_ensembl_version <- ds[[temp_ensembl]]
        }

        temp <- list(
            id              = d,
            annotations     = annotations,
            covar_info      = covar_info,
            datatype        = ds$datatype,
            display_name    = display_name,
            ensembl_version = ds_ensembl_version,
            samples         = ds_synchronized$samples,
            sample_id_field = get_sample_id_field(ds)
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
    datasets <- utils::apropos('^dataset\\.*', ignore.case = TRUE)
    ret <- c()

    for (d in datasets) {
        ds <- get(d)

        annots_field <- NA

        if (ds$datatype == 'mrna') {
            annots_field <- grep("^annots?(\\.|_){1}mrnas?$",
                                 names(ds),
                                 value = TRUE)
        } else if(ds$datatype == 'protein') {
            annots_field <- grep("^annots?(\\.|_){1}proteins?$",
                                 names(ds),
                                 value = TRUE)
        } else if(ds$datatype == 'phos') {
            annots_field <- grep("^annots?(\\.|_){1}phos?$",
                                 names(ds),
                                 value = TRUE)
        } else if(ds$datatype == 'phenotype') {
            annots_field <- grep("^annots?(\\.|_){1}pheno(type)?s?$",
                                 names(ds),
                                 value = TRUE)
        }

        annots_field_samples <- grep("^annots?(\\.|_){1}samples?$",
                                     names(ds),
                                     value = TRUE)

        display_name_field <- grep(
            "^display(\\.|_){1}name$",
            names(ds),
            ignore.case = TRUE,
            value = TRUE
        )

        display_name <- d

        if (length(display_name_field) != 0) {
            display_name <- ds[[display_name_field]]
        }

        temp <- list(id              = d,
                     display_name    = display_name,
                     datatype        = ds$datatype,
                     num_annotations = NROW(ds[[annots_field]]),
                     num_samples     = NROW(ds[[annots_field_samples]]))

        ret <- c(ret, list(temp))
    }

    ret
}

#' Check if id exists and has data in dataset.
#'
#' @param id the id to check
#' @param ds a dataset object
#'
#' @return `TRUE` if id contains data, `FALSE` otherwise
#' @export
id_exists <- function(id, ds = NULL) {
    ret <- list()

    if (valid(ds)) {
        if (is.character(ds)) {
            ds <- get_dataset_by_id(ds)
        }
        ds <- synchronize_dataset(ds)

        if (ds$datatype == 'mrna') {
            all_ids <- ds$annot_mrna$gene_id
        } else if(ds$datatype == 'protein') {
            all_ids <- ds$annot_protein$gene_id
        } else if(ds$datatype == 'phos') {
            all_ids <- ds$annot_phos$gene_id
        } else if(ds$datatype == 'phenotype') {
            all_ids <- ds$annot_phenotype$data_name
        }

        if (any(id == colnames(ds$data)) && (id %in% all_ids)) {
            ret[[d]] <- list(
                dataset_id   = d,
                dataset_name = ds$display_name,
                id           = id
            )
        }
    } else {
        # check all datasets
        datasets <- utils::apropos('^dataset\\.*', ignore.case = TRUE)

        for (d in datasets) {
            ds <- synchronize_dataset(get_dataset_by_id(d))
            all_ids <- NULL

            if (ds$datatype == 'mrna') {
                all_ids <- ds$annot_mrna$gene_id
            } else if(ds$datatype == 'protein') {
                all_ids <- ds$annot_protein$gene_id
            } else if(ds$datatype == 'phos') {
                all_ids <- ds$annot_phos$gene_id
            } else if(ds$datatype == 'phenotype') {
                all_ids <- ds$annot_phenotype$data_name
            }

            if (any(id == colnames(ds$data)) && (id %in% all_ids)) {
                ret[[d]] <- list(
                    dataset_id   = d,
                    dataset_name = ds$display_name,
                    id           = id
                )
            }
        }
    }

    if (invalid(ret)) {
        ret <- NULL
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
    if (any("datatype" == names(ds))) {
        if (startsWith(tolower(ds$datatype), "pheno")) {
            return(TRUE)
        }
    } else {
        message("datatype not found in dataset")
    }

    FALSE
}


