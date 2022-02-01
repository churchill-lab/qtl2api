#' Get the LOD peaks
#'
#' @param ds the dataset object
#' @param intcovar the interactive covariate
#'
#' @return a `data.frame` with the following columns: marker, chr, bp
#' and depending upon dataset$dataType the following columns:
#' mRNA = gene_id, symbol, gene_chrom, middle, lod
#' protein = protein_id, gene_id, symbol, gene_chrom, middle, lod
#' phenotype = data_name, short_name, description, lod
#'
#' @export
get_lod_peaks <- function(ds, intcovar = NULL) {
    # these functions should not be synchronized
    if (!gtools::invalid(ds$is_synchronized)) {
        stop("dataset should not be synchronized")
    }

    peaks <- NULL

    if (is.null(intcovar)) {
        peaks <- ds$lod.peaks$additive
    } else {
        # find the covar and get the name of the lod peaks
        covar_info <- ds$covar.info %>% janitor::clean_names()

        if (any(intcovar == covar_info$sample_column)) {
            n <- covar_info[covar_info$sample_column == intcovar, ]
            peaks <- ds$lod.peaks[[n$lod_peaks]]
        }
    }

    if (gtools::invalid(peaks)) {
        stop(sprintf("No peaks found for intcovar '%s' in lod.peaks", intcovar))
    }

    # this is a little extra work because we are trying to be nice for users
    # who separate with '.' or '_'
    peaks %<>% janitor::clean_names()
    markers_cleaned <- markers %>% janitor::clean_names()

    # convert from Mbp to bp
    if (all(markers_cleaned$pos < 1000)) {
        markers_cleaned$pos <- markers_cleaned$pos * 1000000
    }

    if (tolower(ds$datatype) == "mrna") {
        annots <- ds$annot.mrna %>% janitor::clean_names()

        if (all(annots$start < 1000)) {
            annots$start <- annots$start * 1000000
        }

        if (all(annots$end < 1000)) {
            annots$end <- annots$end * 1000000
        }

        ret <- annots %>%
            dplyr::inner_join(
                peaks,
                by = "gene_id"
            ) %>%
            dplyr::select(
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                gene_chr   = .data$chr,
                start      = .data$start,
                end        = .data$end,
                marker_id  = .data$marker_id,
                lod        = .data$lod
            ) %>%
            dplyr::mutate(
                gene_pos = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                marker_id = .data$marker_id,
                chr       = .data$chr,
                pos       = .data$pos,
                gene_id   = .data$gene_id,
                symbol    = .data$symbol,
                gene_chr  = .data$gene_chr,
                gene_pos  = .data$gene_pos,
                lod       = .data$lod
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
                    by = c("gene_id", "marker_id", "lod")
                )
        }
    } else if (tolower(ds$datatype) == "protein") {
        annots <- ds$annot.protein %>% janitor::clean_names()

        if (all(annots$start < 1000)) {
            annots$start <- annots$start * 1000000
        }

        if (all(annots$end < 1000)) {
            annots$end <- annots$end * 1000000
        }

        ret <- annots %>%
            dplyr::inner_join(
                peaks,
                by = "protein_id"
            ) %>%
            dplyr::select(
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                gene_chr   = .data$chr,
                start      = .data$start,
                end        = .data$end,
                marker_id  = .data$marker_id,
                lod        = .data$lod
            ) %>%
            dplyr::mutate(
                gene_pos = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                marker_id  = .data$marker_id,
                chr        = .data$chr,
                pos        = .data$pos,
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                gene_chr   = .data$gene_chr,
                gene_pos   = .data$gene_pos,
                lod        = .data$lod
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
                    by = c("protein_id", "marker_id", "lod")
                )
        }
    } else if (is_phenotype(ds)) {
        ret <- ds$annot.phenotype %>%
            janitor::clean_names() %>%
            dplyr::inner_join(
                peaks,
                by = "data_name"
            ) %>%
            dplyr::select(
                data_name   = .data$data_name,
                short_name  = .data$short_name,
                description = .data$description,
                marker_id   = .data$marker_id,
                lod         = .data$lod
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                marker_id   = .data$marker_id,
                chr         = .data$chr,
                pos         = .data$pos,
                data_name   = .data$data_name,
                short_name  = .data$short_name,
                description = .data$description,
                lod         = .data$lod
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
        stop(sprintf("dataset has an invalid datatype '%s'", ds$datatype))
    }

    ret
}


#' Get the LOD peaks for additive and all covariates.
#'
#' @param ds the dataset object (not synchronized)
#'
#' @return a data.frame with the following columns: marker, chr, pos
#' and depending upons dataset$dataType the following columns:
#' mRNA = gene_id, symbol, gene_chrom, middle, lod
#' protein = protein_id, gene_id, symbol, gene_chrom, middle, lod
#' phenotype = data_name, short_name, description, lod
#'
#' @export
get_lod_peaks_all <- function(ds) {
    # these functions should not be synchronized
    if (!gtools::invalid(ds$is_synchronized)) {
        stop("dataset should not be synchronized")
    }

    # get the additive LOD peaks
    peaks <- list(additive = get_lod_peaks(ds))

    covar_info <- ds$covar.info %>% janitor::clean_names()

    # get the rest
    for (i in seq(nrow(covar_info))) {
        cov_inf <- covar_info[i, ]

        if (cov_inf$interactive == TRUE) {
            peaks[[cov_inf$sample_column]] <-
                get_lod_peaks(ds, cov_inf$sample_column)
        }
    }

    peaks
}

