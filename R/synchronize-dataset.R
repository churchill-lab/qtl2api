

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
        annot_samples            = ds_synch$samples,
        annots_only_in_data      = ds_synch$annots_only_in_data,
        annots_only_in_annots    = ds_synch$annots_only_in_annots,
        samples_only_in_data     = ds_synch$samples_only_in_data,
        samples_only_in_samples  = ds_synch$samples_only_in_samples,
        covar_info               = covar_info,
        data                     = ds_synch$data,
        datatype                 = datatype,
        display_name             = dataset[[display_name_field]],
        ensembl_version          = ensembl_version,
        sample_id_field          = sample_id_field,
        is_synchronized          = TRUE
    )

    if (tolower(dataset$datatype) == 'mrna') {
        if (all(ds_synch$annots$start < 1000)) {
            ds_synch$annots$start <- as.integer(ds_synch$annots$start * 1000000)
        }

        if (all(ds_synch$annots$end < 1000)) {
            ds_synch$annots$end <- as.integer(ds_synch$annots$end * 1000000)
        }

        ds$annot_mrna <- ds_synch$annots
    } else if (tolower(dataset$datatype) == 'protein') {
        if (all(ds_synch$annots$start < 1000)) {
            ds_synch$annots$start <- as.integer(ds_synch$annots$start * 1000000)
        }

        if (all(ds_synch$annots$end < 1000)) {
            ds_synch$annots$end <- as.integer(ds_synch$annots$end * 1000000)
        }

        ds$annot_protein <- ds_synch$annots
    } else if (tolower(dataset$datatype) == 'protein_uniprot') {
        ds_synch$annots$start <-
            ds_synch$annots$start %>%
            tidyr::replace_na(0)

        if (all(ds_synch$annots$start < 1000)) {
            ds_synch$annots$start <- as.integer(ds_synch$annots$start * 1000000)
        }

        ds_synch$annots$end <-
            ds_synch$annots$end %>%
            tidyr::replace_na(0)

        if (all(ds_synch$annots$end < 1000)) {
            ds_synch$annots$end <- as.integer(ds_synch$annots$end * 1000000)
        }

        ds$annot_protein_uniprot <- ds_synch$annots
    } else if (tolower(dataset$datatype) == 'phos') {
        if (all(ds_synch$annots$start < 1000)) {
            ds_synch$annots$start <- as.integer(ds_synch$annots$start * 1000000)
        }

        if (all(ds_synch$annots$end < 1000)) {
            ds_synch$annots$end <- as.integer(ds_synch$annots$end * 1000000)
        }

        ds$annot_phos <- ds_synch$annots
    } else if (tolower(dataset$datatype) == 'ptm') {
        ds_synch$annots$start <-
            ds_synch$annots$start %>%
            tidyr::replace_na(0)

        if (all(ds_synch$annots$start < 1000)) {
            ds_synch$annots$start <- as.integer(ds_synch$annots$start * 1000000)
        }

        ds_synch$annots$end <-
            ds_synch$annots$end %>%
            tidyr::replace_na(0)

        if (all(ds_synch$annots$end < 1000)) {
            ds_synch$annots$end <- as.integer(ds_synch$annots$end * 1000000)
        }

        ds$annot_ptm <- ds_synch$annots
    } else if (tolower(dataset$datatype) == 'peptide') {
        ds_synch$annots$start <-
            ds_synch$annots$start %>%
            tidyr::replace_na(0)

        if (all(ds_synch$annots$start < 1000)) {
            ds_synch$annots$start <- as.integer(ds_synch$annots$start * 1000000)
        }

        ds_synch$annots$end <-
            ds_synch$annots$end %>%
            tidyr::replace_na(0)

        if (all(ds_synch$annots$end < 1000)) {
            ds_synch$annots$end <- as.integer(ds_synch$annots$end * 1000000)
        }

        ds$annot_peptide <- ds_synch$annots
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
