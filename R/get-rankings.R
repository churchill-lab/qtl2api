#' Calculate the rankings for which gene/protein should be shown.
#'
#' @param ds the dataset object
#' @param chrom the chromosome to filter on
#' @param min_value minimum score value (defaults to 100).
#' @param max_value maximum score value (defaults to 1000).
#'
#' @return A tibble with the following columns: id, ranking
#' @export
get_rankings <- function(ds, chrom = NULL,
                         min_value = 100, max_value = 1000) {
    if (is_phenotype(ds)) {
        stop("phenotype datasets are not supported")
    }

    # different than get_data()
    # raw, norm, rz, log, transformed
    if (is.matrix(ds$data)) {
        data <- ds$data
    } else {
        if (!is.null(ds$data$raw)) {
            data <- ds$data$raw
        } else if (!is.null(ds$data$norm)) {
            data <- ds$data$norm
        } else if (!is.null(ds$data$rz)) {
            data <- ds$data$rz
        } else if (!is.null(ds$data$log)) {
            data <- ds$data$log
        } else if (!is.null(ds$data$transformed)) {
            data <- ds$data$transformed
        }
    }

    # get the mean for each gene/protein/phenotype
    data_mean <- colMeans(data, na.rm = TRUE)

    # grab the min and max values
    min_data_mean <- min(data_mean)
    max_data_mean <- max(data_mean)

    # scale the values between [min_value, max_value]
    tmp <- (data_mean - min_data_mean) / (max_data_mean - min_data_mean)
    tmp <- (max_value - min_value) * tmp + min_value

    if (tolower(ds$datatype) == "mrna") {
        if (!is.null(chrom)) {
            # filter the data to just return the chromosome asked for
            gene_ids <-
                ds$annot.mrna %>%
                janitor::clean_names() %>%
                dplyr::filter(.data$chrom == chrom)

            tmp <- tmp[gene_ids$gene_id]
        }

        ret <- tibble::tibble(
            gene_id = names(tmp),
            ranking = as.integer(tmp)
        )

        # group by gene_id snd than take the gene_id ranking value
        ret <- ret %>%
            dplyr::group_by(.data$gene_id) %>%
            dplyr::summarise(
                ranking = max(.data$ranking)
            )

        return(ret)
    } else {
        annot_protein <- ds$annot.protein %>% janitor::clean_names()

        if (!is.null(chrom)) {
            # filter the data to just return the chromosome asked for
            protein_ids <-
                annot_protein %>%
                dplyr::filter(.data$chrom == chrom)

            tmp <- tmp[protein_ids$protein_id]
        }

        ret <- tibble::tibble(
            id = names(tmp),
            ranking = as.integer(tmp)
        )

        # group by gene_id and than take the gene_id ranking value
        ret <- ret %>%
            dplyr::inner_join(
                annot_protein,
                by = c("id" = "protein_id")
            ) %>%
            dplyr::select(
                protein_id = .data$id,
                gene_id    = .data$gene_id,
                ranking    = .data$ranking
            ) %>%
            dplyr::group_by(.data$gene_id) %>%
            dplyr::summarise(
                ranking = max(.data$ranking)
            )

        return(ret)
    }
}
