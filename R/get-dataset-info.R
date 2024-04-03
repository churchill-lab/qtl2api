#' Get all "dataset.*" information
#'
#' This will return a named list of all datasets and the ensmebl version.
#'
#' @return A named list of all the dataset objects, along with the
#'   ensembl.version.
#' @export
get_dataset_info <- function() {
    datasets <- utils::apropos('^dataset\\.*', ignore.case = TRUE)
    ret <- list()

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

        if (tolower(ds$datatype) == 'mrna') {
            annotations <-
                tibble::tibble(
                    gene_id    = ds_synchronized$annots$gene_id
                )
        } else if(tolower(ds$datatype) == 'protein') {
            annotations <-
                tibble::tibble(
                    protein_id = ds_synchronized$annots$protein_id,
                    gene_id    = ds_synchronized$annots$gene_id
                )
        } else if(tolower(ds$datatype) == 'protein_uniprot') {
            annotations <-
                tibble::tibble(
                    uniprot_id  = ds_synchronized$annots$uniprot_id,
                    protein_id  = ds_synchronized$annots$protein_id,
                    gene_id     = ds_synchronized$annots$gene_id,
                    gene_symbol = ds_synchronized$annots$symbol
                )
        } else if(tolower(ds$datatype) == 'phos') {
            annotations <-
                tibble::tibble(
                    phos_id    = ds_synchronized$annots$phos_id,
                    protein_id = ds_synchronized$annots$protein_id,
                    gene_id    = ds_synchronized$annots$gene_id
                )
        } else if(tolower(ds$datatype) == 'ptm') {
            annotations <-
                tibble::tibble(
                    ptm_id      = ds_synchronized$annots$ptm_id,
                    peptide_id  = ds_synchronized$annots$peptide_id,
                    protein_id  = ds_synchronized$annots$protein_id,
                    gene_id     = ds_synchronized$annots$gene_id,
                    uniprot_id  = ds_synchronized$annots$uniprot_id,
                    gene_symbol = ds_synchronized$annots$symbol
                )

            if('ptm' %in% colnames(ds_synchronized$annots)) {
                annotations['ptm'] = ds_synchronized$annots$ptm
            }
        } else if(tolower(ds$datatype) == 'peptide') {
            annotations <-
                tibble::tibble(
                    peptide_id  = ds_synchronized$annots$peptide_id,
                    protein_id  = ds_synchronized$annots$protein_id,
                    gene_id     = ds_synchronized$annots$gene_id,
                    uniprot_id  = ds_synchronized$annots$uniprot_id,
                    gene_symbol = ds_synchronized$annots$symbol
                )
        } else if(is_phenotype(ds)) {
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
        covar_info_order <- NULL

        if ((length(annots_field) != 0) && (!is.null(ds[[annots_field]]))) {
            covar_info <- ds[[annots_field]] %>% janitor::clean_names()
            covar_info_order <- list()

            for (sample_col in covar_info$sample_column) {
                covar_info_order[[sample_col]] <-
                    levels(ds_synchronized$samples[[sample_col]])
            }
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
            id                       = d,
            annotations              = annotations,
            annots_only_in_data      = ds_synchronized$annots_only_in_data,
            annots_only_in_annots    = ds_synchronized$annots_only_in_annots,
            covar_info               = covar_info,
            covar_info_order         = covar_info_order,
            datatype                 = ds$datatype,
            display_name             = display_name,
            ensembl_version          = ds_ensembl_version,
            samples                  = ds_synchronized$samples,
            sample_id_field          = get_sample_id_field(ds),
            samples_only_in_data     = ds_synchronized$samples_only_in_data,
            samples_only_in_samples  = ds_synchronized$samples_only_in_samples
        )

        ret <- c(ret, list(temp))
    }

    list(datasets        = ret,
         ensembl_version = nvl(ensembl_version, NULL))
}
