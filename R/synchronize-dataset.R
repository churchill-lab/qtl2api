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
    } else if (!is.list(dataset)) {
        stop(paste0("dataset should be a list, but it's a ", class(dataset)))
    }

    if (valid(dataset$is_synchronized)) {
        # already synchronized
        return(dataset)
    }

    # TODO: samples can be top level, but not implemented

    # synchronize the data element
    ds_synch <- synchronize_data(dataset)

    # get the annotatin_info
    annotation_info <- get_annotation_info(dataset)

    # fix the covar_info names
    covar_info <- NULL

    covar_info_name <- grep("^covar?(\\.|_){1}info$",
                            names(dataset),
                            value = TRUE)

    if ((length(covar_info_name) > 0) && (!is.null(dataset[[covar_info_name]]))) {
        covar_info <- dataset[[covar_info_name]]

        if (!is.null(covar_info)) {
            covar_info <- covar_info %>% janitor::clean_names()
        }
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

    display_name_field <- grep(
        "^display(\\.|_){0,1}name$",
        names(dataset),
        ignore.case = TRUE,
        value = TRUE
    )

    datatype <- tolower(dataset$datatype)
    if (startsWith(datatype, "pheno")) {
        datatype <- "phenotype"
    }

    ds <- list(
        annotations                = ds_synch$annotations,
        annotation_info            = annotation_info,
        #annotations_only_in_data   = ds_synch$annotations_only_in_data,
        #annotations_only_in_annots = ds_synch$annotations_only_in_annots,
        covar_info                 = covar_info,
        data                       = ds_synch$data,
        datatype                   = datatype,
        display_name               = dataset[[display_name_field]],
        ensembl_version            = ensembl_version,
        is_synchronized            = TRUE,
        samples                    = ds_synch$samples
        #samples_only_in_data       = ds_synch$samples_only_in_data,
        #samples_only_in_samples    = ds_synch$samples_only_in_samples
    )

    return(ds)
}
