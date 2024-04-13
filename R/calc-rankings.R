#' Calculate the rankings for which gene/protein should be shown.
#'
#' @param dataset The dataset object.
#' @param chrom The chromosome to filter on.
#' @param min_value Minimum score value (defaults to 100).
#' @param max_value Maximum score value (defaults to 1000).
#'
#' @return A `tibble` with the identifiers and rankings.
#'
#' @export
calc_rankings <- function(dataset, chrom = NULL,
                          min_value = 100, max_value = 1000) {
    if (is_phenotype(dataset)) {
        stop("phenotype datasets are not supported")
    }

    # make sure annotations, data, and samples are synchronized
    ds <- synchronize_dataset(dataset)

    # get the mean for each
    data_mean <- colMeans(ds$data, na.rm = TRUE)

    # grab the min and max values
    min_data_mean <- min(data_mean)
    max_data_mean <- max(data_mean)

    # scale the values between [min_value, max_value]
    tmp <- (data_mean - min_data_mean) / (max_data_mean - min_data_mean)
    tmp <- (max_value - min_value) * tmp + min_value

    annotations <- ds$annotations %>% dplyr::filter(.data$gene_id != '')

    if (!is.null(chrom)) {
        annotations <- dplyr::filter(.data$chr == chrom)
    }

    tmp <- tmp[annotations$annotation_id]

    ret <- tibble::tibble(
        id = names(tmp),
        ranking = as.integer(tmp)
    )

    ret <- ret %>%
        dplyr::inner_join(
            ds$annotations,
            by = c("id" = "annotation_id")
        ) %>%
        dplyr::select(
            gene_id    = .data$gene_id,
            ranking    = .data$ranking
        ) %>%
        dplyr::group_by(.data$gene_id) %>%
        dplyr::summarise(
            ranking = max(.data$ranking)
        ) %>%
        dplyr::arrange(dplyr::desc(.data$ranking))

    return(ret)
}
