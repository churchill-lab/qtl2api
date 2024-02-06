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
get_rankings <- function(dataset, chrom = NULL,
                         min_value = 100, max_value = 1000) {
    if (is_phenotype(dataset)) {
        stop("phenotype datasets are not supported")
    }

    # make sure annotations, data, and samples are synchronized
    ds <- synchronize_dataset(dataset)

    # get the mean for each gene/protein/phenotype
    data_mean <- colMeans(ds$data, na.rm = TRUE)

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
                ds$annot_mrna %>%
                dplyr::filter(.data$chr == chrom)

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
    } else if (tolower(ds$datatype) == "protein") {
        if (!is.null(chrom)) {
            # filter the data to just return the chromosome asked for
            protein_ids <-
                ds$annot_protein %>%
                dplyr::filter(.data$chr == chrom)

            tmp <- tmp[protein_ids$protein_id]
        }

        ret <- tibble::tibble(
            id = names(tmp),
            ranking = as.integer(tmp)
        )

        # group by gene_id and than take the gene_id ranking value
        ret <- ret %>%
            dplyr::inner_join(
                ds$annot_protein,
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
    } else if (tolower(ds$datatype) == "protein_uniprot") {
        if (!is.null(chrom)) {
            # filter the data to just return the chromosome asked for
            uniprot_ids <-
                ds$annot_protein_uniprot %>%
                dplyr::filter(.data$chr == chrom)

            tmp <- tmp[uniprot_ids$uniprot_id]
        }

        ret <- tibble::tibble(
            id = names(tmp),
            ranking = as.integer(tmp)
        )

        # group by gene_id and than take the gene_id ranking value
        ret <- ret %>%
            dplyr::inner_join(
                ds$annot_protein,
                by = c("id" = "uniprot_id")
            ) %>%
            dplyr::select(
                uniprot_id = .data$id,
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                ranking    = .data$ranking
            ) %>%
            dplyr::group_by(.data$gene_id) %>%
            dplyr::summarise(
                ranking = max(.data$ranking)
            )

        return(ret)
    } else if (tolower(ds$datatype) == "phos") {
        if (!is.null(chrom)) {
            # filter the data to just return the chromosome asked for
            phos_ids <-
                ds$annot_phos %>%
                dplyr::filter(.data$chr == chrom)

            tmp <- tmp[phos_ids$phos_id]
        }

        ret <- tibble::tibble(
            id = names(tmp),
            ranking = as.integer(tmp)
        )

        # group by gene_id and than take the gene_id ranking value
        ret <- ret %>%
            dplyr::inner_join(
                ds$annot_phos,
                by = c("id" = "phos_id")
            ) %>%
            dplyr::select(
                phos_id    = .data$id,
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                ranking    = .data$ranking
            ) %>%
            dplyr::group_by(.data$gene_id) %>%
            dplyr::summarise(
                ranking = max(.data$ranking)
            )

        return(ret)
    } else if (tolower(ds$datatype) == "ptm") {
        if (!is.null(chrom)) {
            # filter the data to just return the chromosome asked for
            ptm_ids <-
                ds$annot_ptm %>%
                dplyr::filter(.data$chr == chrom)

            tmp <- tmp[ptm_ids$ptm_id]
        }

        ret <- tibble::tibble(
            id = names(tmp),
            ranking = as.integer(tmp)
        )

        # group by uniprot_id and than take the uniprot_id ranking value
        ret <- ret %>%
            dplyr::inner_join(
                ds$annot_ptm,
                by = c("id" = "ptm_id")
            ) %>%
            dplyr::select(
                ptm_id     = .data$id,
                peptide_id = .data$peptide_id,
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                uniprot_id = .data$uniprot_id,
                ranking    = .data$ranking
            ) %>%
            dplyr::group_by(.data$uniprot_id) %>%
            dplyr::summarise(
                ranking = max(.data$ranking)
            )

        return(ret)
    } else if (tolower(ds$datatype) == "peptide") {
        if (!is.null(chrom)) {
            # filter the data to just return the chromosome asked for
            peptide_ids <-
                ds$annot_peptide %>%
                dplyr::filter(.data$chr == chrom)

            tmp <- tmp[peptide_ids$peptide_id]
        }

        ret <- tibble::tibble(
            id = names(tmp),
            ranking = as.integer(tmp)
        )

        # group by uniprot_id and than take the uniprot_id ranking value
        ret <- ret %>%
            dplyr::inner_join(
                ds$annot_ptm,
                by = c("id" = "peptide_id")
            ) %>%
            dplyr::select(
                peptide_id = .data$id,
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                uniprot_id = .data$uniprot_id,
                ranking    = .data$ranking
            ) %>%
            dplyr::group_by(.data$uniprot_id) %>%
            dplyr::summarise(
                ranking = max(.data$ranking)
            )

        return(ret)
    } else {
        stop(paste0(ds$datatype, ' datatype not supported'))
    }
}
