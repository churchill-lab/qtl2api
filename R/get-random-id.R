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
    } else if (tolower(dataset$datatype) == "protein_uniprot") {
        if (invalid(dataset$is_synchronized)) {
            annots_field <- grep("^annots?(\\.|_){1}proteins?(\\.|_){1}uniprots?$",
                                 names(dataset),
                                 value = TRUE)
            annot_ids <-
                dataset[[annots_field]] %>%
                janitor::clean_names()
        } else {
            annot_ids <- dataset$annot_protein_uniprot
        }

        annot_ids <- annot_ids$uniprot_id
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
    } else if (tolower(dataset$datatype) == "ptm") {
        if (invalid(dataset$is_synchronized)) {
            annots_field <- grep("^annots?(\\.|_){1}ptm$",
                                 names(dataset),
                                 value = TRUE)
            annot_ids <-
                dataset[[annots_field]] %>%
                janitor::clean_names()
        } else {
            annot_ids <- dataset$annot_ptm
        }

        annot_ids <- annot_ids$ptm_id
    } else if (tolower(dataset$datatype) == "peptide") {
        if (invalid(dataset$is_synchronized)) {
            annots_field <- grep("^annots?(\\.|_){1}peptide$",
                                 names(dataset),
                                 value = TRUE)
            annot_ids <-
                dataset[[annots_field]] %>%
                janitor::clean_names()
        } else {
            annot_ids <- dataset$annot_peptide
        }

        annot_ids <- annot_ids$peptide_id
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

