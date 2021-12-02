
#' Perform a LOD scan.
#'
#' @param ds the dataset object
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
get_lod_scan <- function(ds, id, intcovar = NULL, cores = 0,
                         filter_threshold = 6.0, filter_peak_drop = Inf,
                         filter_thresholdX = NULL, filter_peak_dropX = NULL,
                         scan1_output = FALSE) {
    # get the data
    data <- get_data(ds)

    # check if id exists
    idx <- which(colnames(data) == id)

    if (gtools::invalid(idx)) {
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
        if (intcovar %not in% ds$covar.info$sample_column) {
            stop(sprintf(
                "intcovar '%s' not found in covar.info in dataset",
                intcovar
            ))
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
        pheno     = data[, id, drop = FALSE],
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

    # construct a 2 dimensional array of data with id, chr, pos, lod as columns
    # we perform a left join here to make sure that the number of elements match
    # also convert from type scan1 to numeric
    lod_scores_mod <-
        dplyr::inner_join(
            tibble::as_tibble(lod_scores, rownames = "marker.id"),
            markers,
            by = "marker.id"
        ) %>%
        dplyr::select(
            id  = .data$marker.id,
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
            markers,
            by = c("chr", "pos")
        ) %>%
        dplyr::select(
            id  = .data$marker.id,
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
#' @param ds the dataset object
#' @param id the unique id in the dataset
#' @param intcovar the interactive covariate
#' @param chrom The chromosome.
#' @param cores number of cores to use (0=ALL)
#'
#' @return A named list with each name a "sample value" and the element is a
#'   tibble with the following columns: id, chr, pos, lod.
#'
#' @importFrom rlang .data
#' @export
get_lod_scan_by_sample <- function(ds, id, intcovar, chrom, cores = 0) {
    # get the data
    data <- get_data(ds)

    # check if id exists
    idx <- which(colnames(data) == id)

    if (gtools::invalid(idx)) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # make sure the chromosome data exists
    if (gtools::invalid(K[[chrom]])) {
        stop(sprintf("Cannot find chromosome '%s' in Kinship matrix", id))
    }

    # extract the markers just for the chromosome
    markers_chrom <- markers %>% dplyr::filter(.data$chr == chrom)

    # make sure nCores is appropriate
    num_cores <- nvl_int(cores, 0)

    if (intcovar %not in% ds$covar.info$sample_column) {
        stop(sprintf(
            "intcovar '%s' not found in covar.info in dataset",
            intcovar
        ))
    }

    # get all the unique values for the interactive.covar and sort them
    if (is.factor(ds$annot.samples[[intcovar]])) {
        covar_unique <- gtools::mixedsort(levels(ds$annot.samples[[intcovar]]))
    } else {
        covar_unique <- gtools::mixedsort(unique(ds$annot.samples[[intcovar]]))
    }

    # convert samples to data.frame because QTL2 relies heavily
    # on rownames and colnames, rownames currently are or will
    # soon be deprecated in tibbles
    samples <- as.data.frame(ds$annot.samples)

    # get the sample id field
    sample_id_field <- get_sample_id_field(ds)

    # set the rownames so scan1 will work
    rownames(samples) <-
        (samples %>% dplyr::select(dplyr::matches(sample_id_field)))[[1]]

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
            ds$annot.samples %>%
            dplyr::filter(!!as.name(intcovar) == u) %>%
            dplyr::select(dplyr::matches(sample_id_field))

        sample_names <- c(sample_names[[1]])

        # get the covar data
        covar <- get_covar_matrix(ds, id)

        # filter by the samples we need
        covar <- ds$covar.matrix[sample_names, , drop = FALSE]

        # exclude covar columns that contain it's name
        # since this is a matrix we cannot use %>% select
        covar <-
            covar[, -which(grepl(intcovar, colnames(covar), ignore.case = T))]

        temp <- qtl2::scan1(
            genoprobs = genoprobs[sample_names, chrom],
            kinship   = K[[chrom]][sample_names, sample_names],
            pheno     = data[sample_names, id, drop = FALSE],
            addcovar  = covar,
            cores     = num_cores,
            reml      = TRUE
        )

        # construct a 2 dimensional array of data with id, chr, pos, lod
        # convert from type scan1 to numeric
        temp <-
            dplyr::inner_join(
                tibble::as_tibble(temp, rownames = "marker.id"),
                markers,
                by = "marker.id"
            ) %>%
            dplyr::select(
                id  = .data$marker.id,
                chr = .data$chr,
                pos = .data$pos,
                lod = id
            ) %>%
            dplyr::mutate_at(c("lod"), as.numeric)

        ret[[toString(u)]] <- temp
    }

    ret
}

