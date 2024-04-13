#' Get all the peaks for an id.
#
#' @param dataset the dataset object
#' @param id the unique id in the dataset
#' @param threshold if set, qtl2::find_peaks is used
#' @param peakdrop if set, qtl2::find_peaks is used
#' @param thresholdX if set, qtl2::find_peaks is used
#' @param peakdropX if set, qtl2::find_peaks is used
#' @param n_cores number of cores to use (0=ALL)
#' @param calc_diff compute the difference between additive and covariate
#'
#' @return a tibble of the peaks.
#'
#' @export
calc_lod_peaks_for_annot <- function(dataset, id,
                                     threshold = 6.0, peakdrop = 2,
                                     thresholdX = 6.0, peakdropX = 2,
                                     n_cores = 0, calc_diff = FALSE) {

    # make sure samples and annotations are available
    ds <- synchronize_dataset(dataset)

    # get the lod scan
    lods_additive <- calc_lod_scores(
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

    if (invalid(peaks_additive)) {
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

            lods_covar <- calc_lod_scores(
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

            if (calc_diff) {
                # DO NOT STORE LOD, STORE LOD DIFF
                temp <- lods_covar$scan1 - lods_additive$scan1
            } else {
                temp <- lods_covar$scan1
            }

            peaks_covar <- qtl2::find_peaks(
                temp,
                map,
                threshold  = threshold,
                peakdrop   = peakdrop,
                thresholdX = thresholdX,
                peakdropX  = peakdropX
            )

            if(valid(peaks_covar)) {
                peaks_covar$scan = inf$sample_column
                peaks_all <- peaks_all %>% dplyr::bind_rows(peaks_covar)
            }
        }
    }

    output <- NULL

    if (nrow(peaks_all) > 0) {
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

        markers_cleaned <- get_markers()

        lod_peaks <-
            dplyr::left_join(
                peaks_all,
                markers_cleaned,
                by = c('chr', 'pos')
            ) %>%
            dplyr::select(
                annotation_id  = .data$lodcolumn,
                marker_id      = .data$marker_id,
                chr            = .data$chr,
                pos            = .data$pos,
                lod            = .data$lod,
                scan           = .data$scan
            ) %>%
            dplyr::mutate_at(
                c("lod"), as.numeric
            ) %>%
            tibble::as_tibble()


        for (i in 1:nrow(lod_peaks)) {
            peak <- lod_peaks[i, ]

            if(invalid(peak$marker_id)) {
                mrk_id <- qtl2::find_marker(
                    map,
                    peak$chr,
                    peak$pos
                )

                mrk <- markers_cleaned %>%
                    dplyr::filter(
                        .data$marker_id == mrk_id
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
