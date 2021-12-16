#' Get the founder coefficients
#'
#' @param dataset the dataset object
#' @param id the unique id in the dataset
#' @param chrom the chromosome
#' @param intcovar the interactive covariate
#' @param blup whether or not to perform BLUP
#' @param center whether or not to center the data
#' @param cores number of cores to use (0=ALL)
#'
#' @return a named `list` with each element being a tibble with the following
#'         columns: id, chr, pos, and A-H
#'
#' @importFrom rlang .data
#' @export
get_founder_coefficients <- function(dataset, id, chrom, intcovar = NULL,
                                     blup = FALSE, center = TRUE, cores = 0) {
    # make sure annotations, data, and samples are synchronized
    ds <- synchronize_dataset(dataset)

    # check if id exists
    if (id %not in% colnames(ds$data)) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # make sure the chromosome data exists
    if (gtools::invalid(K[[chrom]])) {
        stop(sprintf("Cannot find chromosome '%s' in Kinship matrix", id))
    }

    # make sure ncores is appropriate
    num_cores <- nvl_int(cores, 0)

    # get the covar data
    covar <- get_covar_matrix(ds, id)

    # this is a little extra work because we are trying to be nice for users
    # who separate with '.' or '_'
    markers_cleaned <- markers %>% janitor::clean_names()

    ret <- list()

    if (gtools::invalid(intcovar)) {
        if (blup) {
            temp <- qtl2::scan1blup(
                genoprobs = genoprobs[, chrom],
                pheno     = ds$data[, id, drop = FALSE],
                kinship   = K[[chrom]],
                addcovar  = covar,
                cores     = num_cores
            )
        } else {
            temp <- qtl2::scan1coef(
                genoprobs = genoprobs[, chrom],
                pheno     = ds$data[, id, drop = FALSE],
                kinship   = K[[chrom]],
                addcovar  = covar
            )
        }

        if (center) {
            a2h <- LETTERS[1:8]
            temp[, a2h] <- temp[, a2h] - rowMeans(temp[, a2h], na.rm = TRUE)
        }

        ret[["additive"]] <-
            dplyr::inner_join(
                tibble::as_tibble(temp, rownames = "marker_id"),
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                id  = .data$marker_id,
                chr = .data$chr,
                pos = .data$pos,
                LETTERS[1:8]
            ) %>%
            dplyr::mutate_at(dplyr::vars(-c("id", "chr", "pos")), as.numeric)
    } else {
        if (intcovar %not in% ds$covar_info$sample_column) {
            stop(sprintf("intcovar '%s' not found in covar.info", intcovar))
        }

        # get all the unique values for the intcovar and sort them
        if (is.factor(ds$annot_samples[[intcovar]])) {
            covar_unique <-
                gtools::mixedsort(levels(ds$annot_samples[[intcovar]]))
        } else {
            covar_unique <-
                gtools::mixedsort(unique(ds$annot_samples[[intcovar]]))
        }

        # convert samples to data.frame because QTL2 relies heavily
        # on rownames and colnames, rownames currently are or will
        # soon be deprecated in tibbles
        samples <- as.data.frame(ds$annot_samples)

        # set the rownames so scan1 will work
        rownames(samples) <-
            (samples %>% dplyr::select(dplyr::matches(ds$sample_id_field)))[[1]]

        # loop through the unique values for the interactive.covar
        for (u in covar_unique) {
            # sample_names will contain ONLY the samples that match x
            # take all samples
            # filter rows by value, i.e. sex = "F
            # select just need the sample id field column
            sample_names <-
                samples %>%
                dplyr::filter(!!as.name(intcovar) == u) %>%
                dplyr::select(dplyr::matches(ds$sample_id_field))

            sample_names <- c(sample_names[[1]])

            # subset to the intersecting data
            sample_names <- intersect(sample_names, rownames(covar))

            # exclude covar columns that contain it's name
            covar_subset <-
                covar[
                    sample_names,
                    -which(grepl(intcovar, colnames(covar), ignore.case = T))
                ]

            if (blup) {
                temp <- qtl2::scan1blup(
                    genoprobs = genoprobs[sample_names, chrom],
                    pheno     = ds$data[sample_names, id, drop = FALSE],
                    kinship   = K[[chrom]][sample_names, sample_names],
                    addcovar  = covar_subset,
                    cores     = num_cores
                )
            } else {
                temp <- qtl2::scan1coef(
                    genoprobs = genoprobs[sample_names, chrom],
                    pheno     = ds$data[sample_names, id, drop = FALSE],
                    kinship   = K[[chrom]][sample_names, sample_names],
                    addcovar  = covar_subset
                )
            }

            if (center) {
                a2h <- LETTERS[1:8]
                temp[, a2h] <- temp[, a2h] - rowMeans(temp[, a2h], na.rm = TRUE)
            }

            ret[[toString(u)]] <-
                dplyr::inner_join(
                    tibble::as_tibble(temp, rownames = "marker_id"),
                    markers_cleaned,
                    by = "marker_id"
                ) %>%
                dplyr::select(
                    id  = .data$marker_id,
                    chr = .data$chr,
                    pos = .data$pos,
                    LETTERS[1:8]
                ) %>%
                dplyr::mutate_at(
                    dplyr::vars(-c("id", "chr", "pos")),
                    as.numeric
                )
        }
    }

    ret
}
