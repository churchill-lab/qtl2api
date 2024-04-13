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

    annotations <- NULL
    annotations_orig <- NULL

    annotations_name <- grep("^annots?$|annotations?$",
                             names(ds_orig),
                             value = TRUE)

    if (length(annotations_name) == 0) {
        message("ERROR   : annotations not found in dataset")
        return()
    }

    if (!tibble::is_tibble(ds_orig[[annotations_name]])) {
        message("ERROR   : annot_mrna should be a tibble, but found: ", class(ds[[annots_field]]))
    }

    annotations_orig <- ds_orig[[annotations_name]]

    # again allow some play in the column name
    annotations_id_field <- grep(
        "^annots?(\\.|_){0,1}ids?$|^annotations?(\\.|_){0,1}ids?$",
        colnames(annotations_orig),
        value = TRUE,
        ignore.case = TRUE
    )

    if (length(annotations_id_field) == 0) {
        message("ERROR   : annotation_id not found in annotations")
        return()
    }

    if (is_phenotype(dataset)) {
        column_field <- grep(
            "^is(\\.|_){1}id$",
            colnames(annotations_orig),
            value = TRUE,
            ignore.case = TRUE
        )

        if (length(column_field) == 0) {
            message("ERROR   : is_id not found in annotations")
            return()
        }

        column_field <- grep(
            "^data(\\.|_){1}name$",
            colnames(annotations_orig),
            value = TRUE,
            ignore.case = TRUE
        )

        if (length(column_field) == 0) {
            message("ERROR   : data_name not found in annotations")
            return()
        }
    }

    ds <- synchronize_dataset(ds_orig)

    annotations <- ds$annotations

    if (any(duplicated(annotations$annotation_id))) {
        message("ERROR   : There are duplicated annotation_ids in annotations")
    }


    if (is_phenotype(ds)) {
        expected_names <- c('data_name', 'short_name', 'description', 'is_id',
                            'category', 'is_numeric', 'is_date', 'is_factor',
                            'factor_levels', 'is_pheno', 'omit', 'use_covar')

        for (n in expected_names) {
            if (n %not in% names(annotations)) {
                message("ERROR   :", n, "not found in annotations")
            }
        }

        expected_logical <- c('is_id', 'is_pheno', 'is_numeric', 'is_date',
                              'is_factor', 'omit')

        for (n in expected_logical) {
            if(!is.logical(annotations[[n]])) {
                message("ERROR   : annotations$", n, " should be logical, not ", class(annots[[n]]))
            }
        }

        column_field <- grep(
            "^is(\\.|_){1}id$",
            colnames(annotations_orig),
            value = TRUE,
            ignore.case = TRUE
        )

        column_field_data_name <- grep(
            "^data(\\.|_){1}name$",
            colnames(annotations_orig),
            value = TRUE,
            ignore.case = TRUE
        )

        if(is.logical(annotations_orig[[column_field]])) {
            # one and only 1 ID
            theID <- annotations_orig[which(annotations_orig[[column_field]] == TRUE),][[column_field_data_name]]

            if(length(theID) != 1) {
                message("ERROR   : annotations$is_id should have 1 and only 1 row set to TRUE")
            }
        }
    } else {

        expected_names <-  c('gene_id', 'symbol', 'chr', 'start', 'end')

        for (n in expected_names) {
            if (n %not in% names(annotations)) {
                message("ERROR   : ", n, " not found in ", annotations)
            }
        }
    }

    num_annots_orig <- NROW(annotations_orig)
    num_annots_synch <- NROW(annotations)

    if (num_annots_orig != num_annots_synch) {
        cat("WARNING : synchronizing annotations, # samples changed from", num_annots_orig, "to", num_annots_synch, "\n")
    }
}
