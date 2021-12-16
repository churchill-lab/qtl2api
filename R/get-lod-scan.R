
#' Perform a LOD scan.
#'
#' @param dataset the dataset object
#' @param id the unique id in the dataset
#' @param intcovar the interactive covariate
#' @param cores number of cores to use (0=ALL)
#' @param filter_threshold if set, qtl2::find_peaks is used
#' @param filter_peak_drop if set, qtl2::find_peaks is used
#' @param filter_thresholdX if set, qtl2::find_peaks is used
#' @param filter_peak_dropX if set, qtl2::find_peaks is used
#' @param scan1_output if `TRUE`, original `qtl2::scan1` data is included in
#'     return
#'
#' @return a `list` with the following elements:
#'   lod_peaks - list of peaks
#'   lod_scores - tibble with the following columns: id, chr, pos, lod
#'   scan1 - `qtl2::scan1` output
#'
#' @importFrom rlang .data
#' @export
get_lod_scan <- function(dataset, id, intcovar = NULL, cores = 0,
                         filter_threshold = 6.0, filter_peak_drop = Inf,
                         filter_thresholdX = NULL, filter_peak_dropX = NULL,
                         scan1_output = FALSE) {
    # make sure samples and annotations are available
    ds <- synchronize_dataset(dataset)

    # check if id exists
    if (id %not in% colnames(ds$data)) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # make sure num_cores is appropriate
    num_cores <- nvl_int(cores, 0)

    # get the covar data
    covar <- get_covar_matrix(ds, id)

    # set the interactive.covariate, to be used in scan1
    # as scan1(intcovar=interactive_covariate)
    if (gtools::invalid(intcovar)) {
        interactive_covariate <- NULL
    } else {
        if (intcovar %not in% ds$covar_info$sample_column) {
            stop(sprintf("intcovar '%s' not found in covar_info", intcovar))
        }

        # grabbing all the columns from covar (covar.matrix) that
        # match, i.e., "batch" will match "batch2", "BATCH3", etc
        interactive_covariate <-
            covar[, which(grepl(intcovar, colnames(covar), ignore.case = T))]
    }

    # perform the scan using QTL2,
    # - addcovar should always be ALL covars
    # - intcovar should be just the interactive covariate column
    lod_scores <- qtl2::scan1(
        genoprobs = genoprobs,
        kinship   = K,
        pheno     = ds$data[, id, drop = FALSE],
        addcovar  = covar,
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

    markers_cleaned <- markers %>% janitor::clean_names()

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

    ret <- list(
        lod_peaks  = lod_peaks,
        lod_scores = lod_scores_mod
    )

    ret$scan1 <- if (scan1_output) lod_scores else NULL

    ret
}


#' Perform the LOD scan for each "value" of the interactive covariate.
#'
#' @param dataset the dataset object
#' @param id the unique id in the dataset
#' @param chrom The chromosome.
#' @param intcovar the interactive covariate
#' @param cores number of cores to $se (0=ALL)
#'
#' @return A named list with each name a "sample value" and the element is a
#'   tibble with the following columns: id, chr, pos, lod.
#'
#' @importFrom rlang .data
#' @export
get_lod_scan_by_sample <- function(dataset, id, chrom, intcovar, cores = 0) {
    # make sure annotations, data, and samples are synchronized
    ds <- synchronize_dataset(dataset)

    # check if id exists
    if (id %not in% colnames(ds$data)) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # make sure the chromosome data exists
    if (gtools::invalid(K[[chrom]])) {
        stop(sprintf("Cannot find chromosome '%s' in Kinship matrix", chrom))
    }

    # make sure nCores is appropriate
    num_cores <- nvl_int(cores, 0)

    if (intcovar %not in% ds$covar_info$sample_column) {
        stop(sprintf("intcovar '%s' not found in covar_info", intcovar))
    }

    # get all the unique values for the interactive.covar and sort them
    if (is.factor(ds$annot_samples[[intcovar]])) {
        covar_unique <- gtools::mixedsort(levels(ds$annot_samples[[intcovar]]))
    } else {
        covar_unique <- gtools::mixedsort(unique(ds$annot_samples[[intcovar]]))
    }

    # convert samples to data.frame because QTL2 relies heavily
    # on rownames and colnames, rownames currently are or will
    # soon be deprecated in tibbles
    samples <- as.data.frame(ds$annot_samples)

    # set the rownames so scan1 will work
    rownames(samples) <-
        (samples %>% dplyr::select(dplyr::matches(ds$sample_id_field)))[[1]]

    markers_cleaned <- markers %>% janitor::clean_names()

    # ret will be a named list of tibbles with LOD scores
    # each name is a unique sample value
    ret <- list()

    # loop through the unique values for the intcovar
    for (u in covar_unique) {
        # samples.names will contain ONLY the samples that match x
        # take all samples
        # filter rows by value, i.e. sex = "F"
        # select just the sample id field  column
        sample_names <-
            ds$annot_samples %>%
            dplyr::filter(!!as.name(intcovar) == u) %>%
            dplyr::select(dplyr::matches(ds$sample_id_field))

        sample_names <- c(sample_names[[1]])

        # get the covar data
        covar <- get_covar_matrix(ds, id)

        # subset to the intersecting data
        sample_names <- intersect(sample_names, rownames(covar))

        # filter by the samples we need
        covar <- covar[sample_names, , drop = FALSE]

        # exclude covar columns that contain it's name
        # since this is a matrix we cannot use %>% select
        covar <-
            covar[, -which(grepl(intcovar, colnames(covar), ignore.case = T))]

        temp <- qtl2::scan1(
            genoprobs = genoprobs[sample_names, chrom],
            kinship   = K[[chrom]][sample_names, sample_names],
            pheno     = ds$data[sample_names, id, drop = FALSE],
            addcovar  = covar,
            cores     = num_cores,
            reml      = TRUE
        )

        # construct a 2 dimensional array of data with id, chr, pos, lod
        # convert from type scan1 to numeric
        temp <-
            dplyr::inner_join(
                tibble::as_tibble(temp, rownames = "marker_id"),
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                id  = .data$marker_id,
                chr = .data$chr,
                pos = .data$pos,
                lod = id
            ) %>%
            dplyr::mutate_at(c("lod"), as.numeric)

        ret[[toString(u)]] <- temp
    }

    ret
}

