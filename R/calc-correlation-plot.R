#' Get the correlation data for plotting.
#'
#' Perform a correlation analysis on an `id` for a `dataset`.  By default, the
#' correlation will be against the same dataset, but a different dataset can
#' be specified with `dataset_correlate`.  Return x and y values for plotting.
#'
#' @param dataset The dataset object
#' @param id The unique id in the dataset.
#' @param dataset_correlate The dataset to correlate against, dataset if INVALID
#'   or NULL.
#' @param id_correlate The identifier from the correlate dataset.
#' @param intcovar The interactive covariate.
#' @param use_qr `TRUE` to use qr decomposition.
#' @param backfill_NA `TRUE` if we need to impute, we backfill the NAs.
#'
#' @return A named `list` with keys of `data` and `datatypes`.
#'
#' @export
calc_correlation_plot <- function(dataset, id,
                                  dataset_correlate, id_correlate,
                                  intcovar = NULL, use_qr = TRUE,
                                  backfill_NA = TRUE) {
    # make sure samples and annotations are available
    ds <- synchronize_dataset(dataset)

    # get the dataset we are correlating against
    if (is.null(dataset_correlate)) {
        ds_correlate <- ds
    } else {
        ds_correlate <- synchronize_dataset(dataset_correlate)
    }

    # make sure we have the same samples
    samples <- intersect(rownames(ds$data), rownames(ds_correlate$data))

    if (length(samples) == 0) {
        stop("No matching samples to correlate")
    }

    data <- ds$data[samples, ]
    data_correlate <- ds_correlate$data[samples, ]

    # check if id exists
    if (!any(id == colnames(data))) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    } else if (!any(id_correlate == colnames(data_correlate))) {
        stop(sprintf("Cannot find id '%s' in dataset", id_correlate))
    }

    # list of sample names that have been imputed
    samples_imputed <- NULL
    samples_not_imputed <- NULL
    covar_formula <- NULL

    if (is.null(intcovar)) {
        x <- data[, id]
        y <- data_correlate[, id_correlate]
    } else {
        # get the covar information
        covar_information <- get_covar_matrix(ds, id)
        covar_matrix <- covar_information$covar_matrix
        covar_formula <- covar_information$covar_formula

        interactive_covariate <-
            colnames(covar_matrix)[grepl(intcovar, colnames(covar_matrix), ignore.case = T)]

        # determine which samples are imputed
        samples_imputed <- unique(names(which(is.na(data[, id]), arr.ind = TRUE)))
        samples_not_imputed <- setdiff(rownames(data[, id]), samples_imputed)

        # find out where the NA's are in data for the id
        na_index  <- which(is.na(data[, id]), arr.ind = TRUE)

        data <-
            calc_residual_matrix(
                variable_matrix    = data,
                adjust_matrix      = covar_matrix,
                variables_interest = c(id),
                variables_compare  = interactive_covariate,
                use_qr             = use_qr
            )

        # backfill the NA's
        if (backfill_NA) {
            data[na_index] <- NA
        }

        # find out where the NA's are for correlation data
        na_index  <- which(is.na(data_correlate), arr.ind = TRUE)

        data_correlate <-
            calc_residual_matrix(
                variable_matrix    = data_correlate,
                adjust_matrix      = covar_matrix,
                variables_interest = colnames(data_correlate),
                variables_compare  = interactive_covariate,
                use_qr             = use_qr
            )

        # backfill the NA's
        if (backfill_NA) {
            data_correlate[na_index] <- NA
        }

        x <- data[, 1]
        y <- data_correlate[, id_correlate]
    }

    # get the intersecting samples and indices
    samples <- intersect(rownames(data), rownames(data_correlate))
    samples_idx <- which(ds$annot_samples[[ds$sample_id_field]] %in% samples)

    # get the covar factors and their data levels
    sample_info <- list()
    datatypes <- NULL

    if (!is.null(ds$covar_info)) {
        datatypes <- list()

        for (s in ds$covar_info$sample_column) {
            stopifnot(!is.null(ds$annot_samples[[s]]))
            sample_info[[toString(s)]] <- ds$annot_samples[samples_idx, ][[s]]

            if (is.factor(ds$annot_samples[[s]])) {
                datatypes[[toString(s)]] <-
                    gtools::mixedsort(levels(ds$annot_samples[[s]]))
            } else {
                datatypes[[toString(s)]] <-
                    gtools::mixedsort(unique(ds$annot_samples[[s]]))
            }
        }

        correlation_plot_data <-
            tibble::as_tibble(
                data.frame(
                    sample_id = rownames(data)[samples_idx],
                    x         = x[samples_idx],
                    y         = y[samples_idx],
                    sample_info,
                    stringsAsFactors = FALSE
                )) %>%
            dplyr::filter(!is.na(.data$x)) %>%  # don't return the NAs
            dplyr::filter(!is.na(.data$y))      # don't return the NAs
    } else {
        correlation_plot_data <-
            tibble::as_tibble(
                data.frame(
                    sample_id = rownames(data)[samples_idx],
                    x         = x[samples_idx],
                    y         = y[samples_idx],
                    stringsAsFactors = FALSE
                )) %>%
            dplyr::filter(!is.na(.data$x)) %>%  # don't return the NAs
            dplyr::filter(!is.na(.data$y))      # don't return the NAs
    }

    ret <- list(
        datatypes     = datatypes,
        data          = correlation_plot_data
    )

    attr(ret, 'covar_formula') <- covar_formula

    ret
}

