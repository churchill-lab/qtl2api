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
    } else if (tolower(dataset$datatype) == 'protein_uniprot') {
        annots_field <- grep("^annots?(\\.|_){1}proteins?(\\.|_){1}uniprots?$",
                             names(dataset),
                             value = TRUE)

        annots <-
            dataset[[annots_field]] %>%
            janitor::clean_names()

        annot_ids <- annots$uniprot_id
    } else if (tolower(dataset$datatype) == 'phos') {
        annots_field <- grep("^annots?(\\.|_){1}phos?$",
                             names(dataset),
                             value = TRUE)

        annots <-
            dataset[[annots_field]] %>%
            janitor::clean_names()

        annot_ids <- annots$phos_id
    } else if (tolower(dataset$datatype) == 'ptm') {
        annots_field <- grep("^annots?(\\.|_){1}ptm?$",
                             names(dataset),
                             value = TRUE)

        annots <-
            dataset[[annots_field]] %>%
            janitor::clean_names()

        annot_ids <- annots$ptm_id
    } else if (tolower(dataset$datatype) == 'peptide') {
        annots_field <- grep("^annots?(\\.|_){1}peptide?$",
                             names(dataset),
                             value = TRUE)

        annots <-
            dataset[[annots_field]] %>%
            janitor::clean_names()

        annot_ids <- annots$peptide_id
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
    samples_only_in_data <- setdiff(
        rownames(data),
        dataset[[annots_field]][[sample_id_field]]
    )
    samples_only_in_samples <- setdiff(
        dataset[[annots_field]][[sample_id_field]],
        rownames(data)
    )

    if (length(sample_ids) == 0) {
        message("There are no samples in common")
    }

    # get the intersecting annotation ids
    annot_ids <- intersect(
        colnames(data),
        annot_ids
    )
    annots_only_in_data <- setdiff(
        colnames(data),
        annot_ids
    )
    annots_only_in_annots <- setdiff(
        annot_ids,
        colnames(data)
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
    } else if (tolower(dataset$datatype) == 'protein_uniprot') {
        annots <- annots %>% dplyr::filter(.data$uniprot_id %in% annot_ids)
    } else if (tolower(dataset$datatype) == 'phos') {
        annots <- annots %>% dplyr::filter(.data$phos_id %in% annot_ids)
    } else if (tolower(dataset$datatype) == 'ptm') {
        annots <- annots %>% dplyr::filter(.data$ptm_id %in% annot_ids)
    } else if (tolower(dataset$datatype) == 'peptide') {
        annots <- annots %>% dplyr::filter(.data$peptide_id %in% annot_ids)
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
        annots                  = annots,
        annots_only_in_data     = annots_only_in_data,
        annots_only_in_annots   = annots_only_in_annots,
        data                    = data,
        samples                 = samples,
        samples_only_in_data    = samples_only_in_data,
        samples_only_in_samples = samples_only_in_samples,
        sample_id_field         = sample_id_field
    )
}
