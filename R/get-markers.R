#' Get the markers (more for API benefits).  Not for qtl2.
#'
#' @param chrom the chromosome to filter on
#'
#' @return A tibble with the following columns: marker_id, chr, bp
#' @export
get_markers <- function(chrom = NULL) {
    ret <- markers %>%
        janitor::clean_names() %>%
        dplyr::select(
            marker_id = .data$marker_id,
            chr       = .data$chr,
            pos       = .data$pos
        )

    if (all(ret$pos < 1000)) {
        # convert from Mbp to bp
        ret$pos <-  ret$pos * 1000000
    }

    if (!is.null(chrom)) {
        ret %<>% dplyr::filter(.data$chr == chrom)
    }

    ret
}

