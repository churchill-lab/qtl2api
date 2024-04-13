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

    for (d in datasets) {
        ds <- get(d)

        # make sure samples and annotations are available
        ds_synchronized <- synchronize_dataset(ds)

        if(is_phenotype(ds)) {
            annotations <- ds_synchronized$annotations
        } else {
            anotation_columns <-
                ds_synchronized$annotation_info %>%
                dplyr::arrange(.data$order)

            anotation_columns <- c('annotation_id', anotation_columns$column)

            annotations <-
                ds_synchronized$annotations %>%
                dplyr::select(
                    anotation_columns,
                    symbol  = .data$symbol,
                    chr     = .data$chr,
                    start   = .data$start,
                    end     = .data$end
                )
        }

        covar_info_order <- list()

        for (sample_col in ds_synchronized$covar_info$sample_column) {
            covar_info_order[[sample_col]] <-
                levels(ds_synchronized$samples[[sample_col]])
        }

        temp <- list(
            id                       = d,
            annotations              = annotations,
            annotation_info          = ds_synchronized$annotation_info,
            #annots_only_in_data      = ds_synchronized$annots_only_in_data,
            #annots_only_in_annots    = ds_synchronized$annots_only_in_annots,
            covar_info               = ds_synchronized$covar_info,
            covar_info_order         = covar_info_order,
            datatype                 = ds_synchronized$datatype,
            display_name             = ds_synchronized$display_name,
            ensembl_version          = ds_synchronized$ensembl_version,
            samples                  = ds_synchronized$samples
            #samples_only_in_data     = ds_synchronized$samples_only_in_data,
            #samples_only_in_samples  = ds_synchronized$samples_only_in_samples
        )

        ret <- c(ret, list(temp))
    }

    ensembl_version <-
        utils::apropos("^ensembl(\\.|_){1}version$", ignore.case = TRUE)

    if (valid(ensembl_version)) {
        ensembl_version <- get(ensembl_version)
    } else {
        ensembl_version <- NULL
    }

    list(
        datasets        = ret,
        ensembl_version = ensembl_version
    )
}
