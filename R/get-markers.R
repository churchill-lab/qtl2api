#' Get the markers (more for API benefits).  Not for qtl2.
#'
#' @param chrom The chromosome to filter on.
#'
#' @return A tibble with the following columns: marker_id, chr, bp
#' @export
get_markers <- function(chrom = NULL) {
    column_names <- colnames(markers)

    marker_id_field <- grep(
        "^markers?(\\.|_){0,1}ids?$",
        column_names,
        value = TRUE,
        ignore.case = TRUE
    )

    chr_field <- grep(
        "^chrs?$|^chroms?$|^chromosomes?$",
        column_names,
        value = TRUE,
        ignore.case = TRUE
    )

    pos_field <- grep(
        "^pos$|positions?$",
        column_names,
        value = TRUE,
        ignore.case = TRUE
    )

    # make sure we pass back a tibble with marker_id, chr, pos
    if (marker_id_field != "marker_id") {
        markers$marker_id <- markers[[marker_id_field]]
    }

    if (chr_field != "chr") {
        markers$chr <- markers[[chr_field]]
    }

    if (pos_field != "pos") {
        markers$pos <- markers[[pos_field]]
    }

    ret <- markers %>%
        dplyr::filter(!is.na(.data$pos)) %>%
        dplyr::select(
            marker_id = .data$marker_id,
            chr       = .data$chr,
            pos       = .data$pos,
        )

    if (all(ret$pos < 1000)) {
        # convert from Mbp to bp
        ret$pos <-  as.integer(ret$pos * 1000000)
    }

    if (!is.null(chrom)) {
        ret %<>% dplyr::filter(.data$chr == chrom)
    }

    ret
}
