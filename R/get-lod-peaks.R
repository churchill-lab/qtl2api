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



#' Get all the peaks for an id.
#
#' @param dataset the dataset object
#' @param id the unique id in the dataset
#' @param threshold if set, qtl2::find_peaks is used
#' @param peakdrop if set, qtl2::find_peaks is used
#' @param thresholdX if set, qtl2::find_peaks is used
#' @param peakdropX if set, qtl2::find_peaks is used
#' @param n_cores number of cores to use (0=ALL)
#'
#' @return a tibble of the peaks.
#'
#' @export
get_lod_peaks_for_annot <- function(dataset, id,
                                    threshold = 6.0, peakdrop = 2,
                                    thresholdX = 6.0, peakdropX = 2,
                                    n_cores = 0) {

    # make sure samples and annotations are available
    ds <- synchronize_dataset(dataset)

    # convert the numeric index into a string id
    if (class(id) ==  "numeric") {
        id <- colnames(ds$data)[id]
    }

    # get the lod scan
    lods_additive <- qtl2api::get_lod_scan(
        ds, id,
        intcovar          = NULL,
        cores             = n_cores,
        filter_threshold  = threshold,
        filter_peak_drop  = peakdrop,
        filter_thresholdX = thresholdX,
        filter_peak_dropX = peakdropX,
        scan1_output      = TRUE
    )

    lods_additive$lod_scores$scan <- 'additive'

    # find the peaks
    peaks_additive <- qtl2::find_peaks(
        lods_additive$scan1,
        map,
        threshold  = threshold,
        peakdrop   = peakdrop,
        thresholdX = thresholdX,
        peakdropX  = peakdropX
    )

    if (gtools::invalid(peaks_additive)) {
        peaks_additive <- tibble::tibble(
            lodindex  = numeric(),
            lodcolumn = character(),
            chr       = character(),
            pos       = numeric(),
            lod       = numeric(),
            scan      = character()
        )
    } else {
        peaks_additive$scan = 'additive'
    }

    peaks_all <- peaks_additive

    # loop through interactive covariates
    for (i in 1:nrow(ds$covar_info)) {
        inf <- ds$covar_info[i, ]
        if (inf$interactive) {

            lods_covar <- qtl2api::get_lod_scan(
                ds, id,
                intcovar          = inf$sample_column,
                cores             = n_cores,
                filter_threshold  = threshold,
                filter_peak_drop  = peakdrop,
                filter_thresholdX = thresholdX,
                filter_peak_dropX = peakdropX,
                scan1_output      = TRUE
            )

            lods_covar$lod_scores$scan <- inf$sample_column

            # DO NOT STORE LOD, SCORE LOD DIFF
            temp <- lods_covar$scan1 - lods_additive$scan1

            peaks_covar <- qtl2::find_peaks(
                temp,
                map,
                threshold  = threshold,
                peakdrop   = peakdrop,
                thresholdX = thresholdX,
                peakdropX  = peakdropX
            )

            if(!gtools::invalid(peaks_covar)) {
                peaks_covar$scan = inf$sample_column
                peaks_all <- peaks_all %>% dplyr::bind_rows(peaks_covar)
            }
        }
    }

    output <- NULL

    if (nrow(peaks_all) > 0) {
        lod_peaks <-
            dplyr::left_join(
                peaks_all,
                markers,
                by = c('chr', 'pos')
            ) %>%
            dplyr::select(
                annot_id  = .data$lodcolumn,
                marker_id = .data$marker.id,
                chr       = .data$chr,
                pos       = .data$pos,
                lod       = .data$lod,
                scan      = .data$scan
            ) %>%
            dplyr::mutate_at(
                c("lod"), as.numeric
            ) %>%
            tibble::as_tibble()

        output <- tibble::tibble(
                annot_id  = character(),
                marker_id = character(),
                chr       = character(),
                pos       = numeric(),
                lod       = numeric(),
                scan      = character(),
                A         = numeric(),
                B         = numeric(),
                C         = numeric(),
                D         = numeric(),
                E         = numeric(),
                F         = numeric(),
                G         = numeric(),
                H         = numeric()
        )

        for (i in 1:nrow(lod_peaks)) {
            peak <- lod_peaks[i, ]

            if(gtools::invalid(peak$marker_id)) {
                mrk_id <- qtl2::find_marker(
                    map,
                    peak$chr,
                    peak$pos
                )

                mrk <- markers %>%
                    dplyr::filter(
                        marker.id == mrk_id
                    )

                peak$marker_id <- mrk_id
                peak$pos <- mrk$pos
                lod_peaks[i, ] <- peak
            }

            # get the allele effects only for additive
            if(peak$scan == 'additive') {
                temp_probs <- list()
                temp_probs[[peak$chr]] <-
                    genoprobs[[peak$chr]][,,peak$marker_id, drop = FALSE]

                af <- qtl2::scan1blup(
                    genoprobs = temp_probs,
                    pheno     = ds$data[, id, drop = FALSE],
                    kinship   = K[[peak$chr]]
                )

                output <- output %>% dplyr::bind_rows(
                    lod_peaks[i, ] %>%
                         dplyr::bind_cols(tibble::as_tibble(t(af[, LETTERS[1:8]])))
                )
            } else {
                output <- output %>% dplyr::bind_rows(
                    lod_peaks[i, ] %>%
                        dplyr::bind_cols(tibble::tibble(A=NA,B=NA,C=NA,D=NA,E=NA,F=NA,G=NA,H=NA))
                )
            }
        }
    }

    output
}


