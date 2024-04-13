#' Get the LOD peaks for a dataset.
#'
#' Used internally to the package.
#'
#' By default, this will get all LOD peaks stored in the `dataset$lod_peaks`
#' `list`.  If `intcovar` is specified, only the peaks for that are returned.
#'
#' @param ds The dataset object to get the LOD peaks for.
#' @param intcovar NULL for 'additive' or a value representing the interactive
#' covariate.
#' @param full_annotations `TRUE` to send back full annotations.
#'
#' @return A a `tibble` with the following columns depending upon
#' `dataset$datatype`.
#'
#' If the allele effects are stored, values `A-H` will also be included.
#' @export
get_lod_peaks_by_covar <- function(ds, intcovar = NULL,
                                   full_annotations = FALSE) {
    # these functions should not be synchronized
    if (valid(ds$is_synchronized)) {
        stop("dataset should not be synchronized")
    }

    annots_field_peaks <- grep("^lod(\\.|_){1}peaks?$",
                               names(ds),
                               value = TRUE)

    if (length(annots_field_peaks) == 0) {
        return(NULL)
    } else if (is.null(ds[[annots_field_peaks]])) {
        return(NULL)
    }

    ds_synch <- synchronize_dataset(ds)

    peaks <- NULL

    if ((is.null(intcovar)) || (intcovar == 'additive')) {
        peaks <- ds[[annots_field_peaks]]$additive
    } else {
        # find the covar and get the name of the lod peaks
        if (any(intcovar == ds_synch$covar_info$sample_column)) {
            n <- ds_synch$covar_info[ds_synch$covar_info$sample_column == intcovar, ]
            peaks <- ds[[annots_field_peaks]][[n$lod_peaks]]
        }
    }

    if (invalid(peaks)) {
        stop(sprintf("No peaks found for intcovar '%s' in lod.peaks", intcovar))
    }

    # this is a little extra work because we are trying to be nice for users
    # who separate with '.' or '_'
    peaks %<>% janitor::clean_names()

    markers_cleaned <- get_markers()

    if (is_phenotype(ds)) {
        ret <- ds_synch$annotations %>%
            dplyr::inner_join(
                peaks,
                by = "annotation_id"
            ) %>%
            dplyr::select(
                annotation_id = .data$annotation_id,
                data_name     = .data$data_name,
                short_name    = .data$short_name,
                description   = .data$description,
                marker_id     = .data$marker_id,
                lod           = .data$lod
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                annotation_id = .data$annotation_id,
                data_name     = .data$data_name,
                short_name    = .data$short_name,
                description   = .data$description,
                marker_id     = .data$marker_id,
                chr           = .data$chr,
                pos           = .data$pos,
                lod           = .data$lod
            ) %>%
            dplyr::arrange(
                .data$chr,
                .data$pos
            )

        # now add A-H for additive if they exist
        if (all(tolower(LETTERS[1:8]) %in% colnames(peaks))) {
            ret <- ret %>%
                dplyr::inner_join(
                    peaks,
                    by = c("data_name", "marker_id", "lod")
                )
        }
    } else {
        vars_to_select <- c("annotation_id")

        if (full_annotations) {
            vars_to_select <- c("annotation_id")
            if (valid(ds_synch$annotation_info)) {
                vars_to_select <- c(vars_to_select, ds_synch$annotation_info$column)
            }
        }

        ret <- ds_synch$annotations %>%
            dplyr::inner_join(
                peaks,
                by = "annotation_id"
            ) %>%
            dplyr::select(
                vars_to_select,
                symbol     = .data$symbol,
                gene_chr   = .data$chr,
                start      = .data$start,
                end        = .data$end,
                marker_id  = .data$marker_id,
                lod        = .data$lod
            ) %>%
            dplyr::mutate(
                gene_pos = as.integer(round((.data$start + .data$end) / 2))
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                vars_to_select,
                lod        = .data$lod,
                marker_id  = .data$marker_id,
                chr        = .data$chr,
                pos        = .data$pos,
                gene_chr   = .data$gene_chr,
                gene_pos   = .data$gene_pos
            ) %>%
            dplyr::arrange(
                .data$chr,
                .data$pos
            )

        # now add A-H for additive if they exist
        if (all(tolower(LETTERS[1:8]) %in% colnames(peaks))) {
            ret <- ret %>%
                dplyr::inner_join(
                    peaks,
                    by = c("annotation_id", "marker_id", "lod")
                )
        }

    }

    ret
}


