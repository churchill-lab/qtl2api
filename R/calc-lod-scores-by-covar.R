
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
get_lod_scores_by_covar <- function(dataset, id, chrom, intcovar, cores = 0) {
    # make sure annotations, data, and samples are synchronized
    ds <- synchronize_dataset(dataset)

    # check if id exists
    if (!any(id == colnames(ds$data))) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # make sure the chromosome data exists
    if (!any(chrom == names(K))) {
        stop(sprintf("Cannot find chromosome '%s' in Kinship matrix", chrom))
    }

    # make sure nCores is appropriate
    num_cores <- nvl_int(cores, 0)

    if (!any(intcovar == ds$covar_info$sample_column)) {
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

    markers_cleaned <-
        markers %>%
        dplyr::filter(!is.na(.data$pos)) %>%
        janitor::clean_names()

    # ret will be a named list of tibbles with LOD scores
    # each name is a unique sample value
    ret <- list()

    # loop through the unique values for the intcovar
    for (u in covar_unique) {
        # samples.names will contain ONLY the samples that match x
        # take all samples
        # filter rows by value, i.e. sex = "F"
        # select just the sample id field column
        sample_names <-
            ds$annot_samples %>%
            dplyr::filter(!!as.name(intcovar) == u) %>%
            dplyr::select(dplyr::matches(ds$sample_id_field))

        sample_names <- c(sample_names[[1]])

        # get the covar data
        covar_information <- get_covar_matrix(ds, id)
        covar_matrix <- covar_information$covar_matrix

        # subset to the intersecting data
        sample_names <- intersect(sample_names, rownames(covar_matrix))
        sample_names <- intersect(sample_names, rownames(K[[chrom]]))

        # filter by the samples we need
        covar_matrix <- covar_matrix[sample_names, , drop = FALSE]

        # exclude covar columns that contain it's name
        # since this is a matrix we cannot use %>% select
        covar_matrix <-
            covar_matrix[, -which(grepl(intcovar, colnames(covar_matrix), ignore.case = T))]

        temp <- qtl2::scan1(
            genoprobs = genoprobs[sample_names, chrom],
            kinship   = K[[chrom]][sample_names, sample_names],
            pheno     = ds$data[sample_names, id, drop = FALSE],
            addcovar  = covar_matrix,
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
            dplyr::mutate_at(
                c("lod"),
                as.numeric
            )

        if (all(temp$pos < 1000)) {
            temp$pos <- temp$pos * 1000000
        }

        attr(temp, 'covar_formula') <- covar_information$covar_formula
        attr(temp, 'samples') <- sample_names

        ret[[toString(u)]] <- temp
    }

    ret
}

