#' Validate the dataset to make sure qtl2api can use it.
#'
#' @param dataset the dataset_id as a string identifier
#' @param by_api TRUE to perform scans on random elements in the data
#'
#' @export
validate_dataset <- function(dataset, by_api = FALSE) {

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

    if (by_api) {
        validate_dataset_by_api(ds_orig)
    }
}
