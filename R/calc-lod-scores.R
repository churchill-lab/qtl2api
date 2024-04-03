#' Perform a LOD scan.
#'
#' @param dataset The dataset object.
#' @param id The unique id in the dataset.
#' @param intcovar The interactive covariate.
#' @param cores Number of cores to use (0 = ALL).
#' @param filter_threshold If set, qtl2::find_peaks is used.
#' @param filter_peak_drop If set, qtl2::find_peaks is used.
#' @param filter_thresholdX If set, qtl2::find_peaks is used.
#' @param filter_peak_dropX If set, qtl2::find_peaks is used.
#' @param scan1_output If `TRUE`, original `qtl2::scan1` data is included.
#'
#' @return a `list` with the following elements:
#' \itemize{
#'   \item lod_peaks - `tibble` of LOD peaks
#'   \item lod_scores - `tibble` with the following columns: id, chr, pos, lod
#'   \item scan1 - `qtl2::scan1` output
#' }
#'
#' @importFrom rlang .data
#' @export
calc_lod_scores <- function(dataset, id, intcovar = NULL, cores = 0,
                            filter_threshold = 6.0, filter_peak_drop = Inf,
                            filter_thresholdX = NULL, filter_peak_dropX = NULL,
                            scan1_output = FALSE) {
    # make sure samples and annotations are available
    ds <- synchronize_dataset(dataset)

    # check if id exists
    if (!any(id == colnames(ds$data))) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # make sure num_cores is appropriate
    num_cores <- nvl_int(cores, 0)

    # get the covar information
    covar_information <- get_covar_matrix(ds, id)
    covar_matrix <- covar_information$covar_matrix
    covar_formula <- covar_information$covar_formula

    # set the interactive.covariate, to be used in scan1
    interactive_covariate <- NULL

    if (!is.null(intcovar)) {
        if (!any(intcovar == ds$covar_info$sample_column)) {
            stop(sprintf("intcovar '%s' not found in covar_info", intcovar))
        }

        if (is.null(covar_matrix)) {
            stop(sprintf("no covar_matrix, but intcovar '%s' specified", intcovar))
        }

        # grabbing all the columns from covar (covar.matrix) that
        # match, i.e., "batch" will match "batch2", "BATCH3", etc
        interactive_covariate <-
            covar_matrix[, which(grepl(intcovar, colnames(covar_matrix), ignore.case = T))]
    }

    # perform the scan using QTL2,
    # - addcovar should always be ALL covars
    # - intcovar should be just the interactive covariate column
    lod_scores <- qtl2::scan1(
        genoprobs = genoprobs,
        kinship   = K,
        pheno     = ds$data[, id, drop = FALSE],
        addcovar  = covar_matrix,
        intcovar  = interactive_covariate,
        cores     = num_cores,
        reml      = TRUE
    )

    # utilize qtl2::find_peaks
    lod_peaks <- qtl2::find_peaks(
        lod_scores,
        map,
        threshold  = filter_threshold,
        peakdrop   = filter_peak_drop,
        thresholdX = filter_thresholdX,
        peakdropX  = filter_peak_dropX
    )

    markers_cleaned <-
        markers %>%
        dplyr::filter(!is.na(.data$pos)) %>%
        janitor::clean_names()

    # construct a 2 dimensional array of data with id, chr, pos, lod as columns
    # we perform a left join here to make sure that the number of elements match
    # also convert from type scan1 to numeric
    lod_scores_mod <-
        dplyr::inner_join(
            tibble::as_tibble(lod_scores, rownames = "marker_id"),
            markers_cleaned,
            by = "marker_id"
        ) %>%
        dplyr::select(
            id  = .data$marker_id,
            chr = .data$chr,
            pos = .data$pos,
            lod = id
        ) %>%
        dplyr::mutate_at(
            c("lod"),
            as.numeric
        )

    lod_peaks <-
        dplyr::inner_join(
            lod_peaks,
            markers_cleaned,
            by = c("chr", "pos")
        ) %>%
        dplyr::select(
            id  = .data$marker_id,
            chr = .data$chr,
            pos = .data$pos,
            lod = .data$lod
        ) %>%
        dplyr::mutate_at(
            c("lod"),
            as.numeric
        ) %>%
        tibble::as_tibble()

    # convert to bp
    if (all(lod_scores_mod$pos < 1000)) {
        lod_scores_mod$pos <- as.integer(lod_scores_mod$pos * 1000000)
    }

    if (all(lod_peaks$pos < 1000)) {
        lod_peaks$pos <- as.integer(lod_peaks$pos * 1000000)
    }

    ret <- list(
        lod_peaks     = lod_peaks,
        lod_scores    = lod_scores_mod
    )

    ret$scan1 <- if (scan1_output) lod_scores else NULL

    attr(ret, 'covar_formula') <- covar_formula

    ret
}
