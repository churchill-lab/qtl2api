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

    annots <- NULL
    annots_orig <- NULL

    # check orginal
    if (tolower(ds_orig$datatype) == "mrna") {
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

        annots_orig <- ds_orig[[annots_field]]
    } else if (tolower(ds_orig$datatype) == "protein") {
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

        annots_orig <- ds_orig[[annots_field]]
    } else if (tolower(ds_orig$datatype) == "protein_uniprot") {
        annots_field <- grep("^annots?(\\.|_){1}proteins?(\\.|_){1}uniprots?$",
                             names(ds_orig),
                             value = TRUE)

        if (length(annots_field) == 0) {
            message("ERROR   : annot_protein_uniprot not found in dataset.")
            return()
        }

        if (!tibble::is_tibble(ds_orig[[annots_field]])) {
            message("ERROR   : annot_protein_uniprot should be a tibble, but found: ", class(ds[[annots_field]]))
        }

        annots_orig <- ds_orig[[annots_field]]
    } else if (tolower(ds_orig$datatype) == "phos") {
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

        annots_orig <- ds_orig[[annots_field]]
    } else if (tolower(ds_orig$datatype) == "ptm") {
        annots_field <- grep("^annots?(\\.|_){1}ptm$",
                             names(ds_orig),
                             value = TRUE)

        if (length(annots_field) == 0) {
            message("ERROR   : annot_ptm not found in dataset.")
            return()
        }

        if (!tibble::is_tibble(ds_orig[[annots_field]])) {
            message("ERROR   : annot_ptm should be a tibble, but found: ", class(ds[[annots_field]]))
        }

        annots_orig <- ds_orig[[annots_field]]
    } else if (tolower(ds_orig$datatype) == "peptide") {
        annots_field <- grep("^annots?(\\.|_){1}peptide$",
                             names(ds_orig),
                             value = TRUE)

        if (length(annots_field) == 0) {
            message("ERROR   : annot_peptide not found in dataset.")
            return()
        }

        if (!tibble::is_tibble(ds_orig[[annots_field]])) {
            message("ERROR   : annot_peptide should be a tibble, but found: ", class(ds[[annots_field]]))
        }

        annots_orig <- ds_orig[[annots_field]]
    } else if (is_phenotype(ds_orig)) {
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

        annots_orig <- ds_orig[[annots_field]]

        column_field <- grep(
            "^is(\\.|_){1}id$",
            colnames(annots_orig),
            value = TRUE,
            ignore.case = TRUE
        )

        if (length(column_field) == 0) {
            message("ERROR   : is_id not found in annots_phenotype")
            return()
        }

        column_field <- grep(
            "^data(\\.|_){1}name$",
            colnames(annots_orig),
            value = TRUE,
            ignore.case = TRUE
        )

        if (length(column_field) == 0) {
            message("ERROR   : data_name not found in annots_phenotype")
            return()
        }
    }

    ds <- synchronize_dataset(ds_orig)

    annots <- NULL

    if (tolower(ds$datatype) == "mrna") {
        annots <- ds$annot_mrna

        if (any(duplicated(annots$gene_id))) {
            message("ERROR   : There are duplicated gene identifiers in annot_mrna")
        }
    } else if (tolower(ds$datatype) == "protein") {
        annots <- ds$annot_protein

        if (any(duplicated(annots$protein_id))) {
            message("ERROR   : There are duplicated protein identifiers in annot_protein")
        }
    } else if (tolower(ds$datatype) == "protein_uniprot") {
        annots <- ds$annot_protein_uniprot

        if (any(duplicated(annots$uniprot_id))) {
            message("ERROR   : There are duplicated uniprot identifiers in annot_protein_uniprot")
        }
    } else if (tolower(ds$datatype) == "phos") {
        annots <- ds$annot_phos

        if (any(duplicated(annots$phos_id))) {
            message("ERROR   : There are duplicated phos identifiers in annot_phos")
        }
    } else if (tolower(ds$datatype) == "ptm") {
        annots <- ds$annot_ptm

        if (any(duplicated(annots$ptm_id))) {
            message("ERROR   : There are duplicated ptm identifiers in annot_ptm")
        }
    } else if (tolower(ds$datatype) == "peptide") {
        annots <- ds$annot_peptide

        if (any(duplicated(annots$peptide_id))) {
            message("ERROR   : There are duplicated peptide identifiers in annot_peptide")
        }
    } else if (is_phenotype(ds)) {
        annots <- ds$annot_phenotype

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
        } else if (tolower(ds$datatype) == "protein_uniprot") {
            annot_name <- 'annot_protein_uniprot'

            if ('uniprot_id' %not in% names(annots)) {
                message("ERROR   : uniprot_id not found in annot_protein_uniprot")
            }

            if ('protein_id' %not in% names(annots)) {
                message("ERROR   : protein_id not found in annot_protein_uniprot")
            }
        } else if (tolower(ds$datatype) == "phos") {
            annot_name <- 'annot_phos'

            if ('protein_id' %not in% names(annots)) {
                message("ERROR   : protein_id not found in annot_phos")
            }

            if ('phos_id' %not in% names(annots)) {
                message("ERROR   : phos_id not found in annot_phos")
            }
        } else if (tolower(ds$datatype) == "ptm") {
            annot_name <- 'annot_ptm'

            if ('ptm_id' %not in% names(annots)) {
                message("ERROR   : ptm_id not found in annot_ptm")
            }
        } else if (tolower(ds$datatype) == "peptide") {
            annot_name <- 'annot_peptide'

            if ('peptide_id' %not in% names(annots)) {
                message("ERROR   : peptide_id not found in annot_peptide")
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
