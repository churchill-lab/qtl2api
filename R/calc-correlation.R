#' Get the correlation.
#'
#' Perform a correlation analysis on an `id` for a `dataset`.  By default, the
#' correlation will be against the same dataset, but a different dataset can
#' be specified with `dataset_correlate`.
#'
#' @param dataset The dataset object
#' @param id The unique id in the dataset.
#' @param dataset_correlate The dataset to correlate against, dataset if INVALID
#'   or NULL.
#' @param intcovar The interactive covariate.
#' @param use_qr `TRUE` to use qr decomposition.
#' @param backfill_NA `TRUE` if we need to impute, we backfill the NAs.
#' @param full_annotations `TRUE` to send back full annotations.
#'
#' @return A `tibble` with the correlation and annotations.
#'
#' @export
calc_correlation <- function(dataset, id, dataset_correlate = NULL,
                             intcovar = NULL, use_qr = TRUE,
                             backfill_NA = TRUE, full_annotations = FALSE) {
    # make sure annotations, data, and samples are synchronized
    ds <- synchronize_dataset(dataset)

    # get the dataset we are correlating against
    if (is.null(dataset_correlate)) {
        ds_correlate <- ds
    } else {
        ds_correlate <- synchronize_dataset(dataset_correlate)
    }

    # check if id exists
    idx <- which(id == colnames(ds$data))
    if (!any(idx)) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # make sure we have the same samples between datasets
    samples <- intersect(rownames(ds$data), rownames(ds_correlate$data))

    if (length(samples) == 0) {
        stop("No matching samples for correlation")
    }

    data <- ds$data[samples, ]
    data_correlate <- ds_correlate$data[samples, ]

    # list of sample names that have been imputed
    samples_imputed <- NULL
    samples_not_imputed <- NULL
    covar_formula <- NULL

    if (is.null(intcovar)) {
        pcor <- stats::cor(data[, id], data_correlate, use = "pair")
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

        id_residual_matrix <-
            calc_residual_matrix(
                variable_matrix    = data,
                adjust_matrix      = covar_matrix,
                variables_interest = c(id),
                variables_compare  = interactive_covariate,
                use_qr             = use_qr
            )

        # backfill the NA's
        if (backfill_NA) {
            id_residual_matrix[na_index] <- NA
        }

        # find out where the NA's are for correlation data
        na_index  <- which(is.na(data_correlate), arr.ind = TRUE)

        residual_matrix <-
            calc_residual_matrix(
                variable_matrix    = data_correlate,
                adjust_matrix      = covar_matrix,
                variables_interest = colnames(data_correlate),
                variables_compare  = interactive_covariate,
                use_qr             = use_qr
            )

        # backfill the NA's
        if (backfill_NA) {
            residual_matrix[na_index] <- NA
        }

        pcor <- stats::cor(
            id_residual_matrix,
            residual_matrix,
            use = "pair"
        )
    }

    # reorder in decreasing order
    pcor <- pcor[1, order(abs(pcor), decreasing = TRUE)]

    correlations <- tibble::tibble(
        annotation_id = names(pcor),
        cor           = pcor
    )

    # attach some annotations to the correlation data
    if (is_phenotype(ds_correlate)) {
        correlations <-
            correlations %>%
            dplyr::inner_join(
                ds_correlate$annotations,
                by = c("annotation_id")
            ) %>%
            dplyr::select(
                snnotation_id = .data$annotation_id,
                short_name    = .data$short_name,
                cor           = .data$cor
            )
    } else {
        vars_to_select <- c("annotation_id", "cor")

        if (full_annotations) {
            vars_to_select <- c("annotation_id")
            if (valid(ds_correlate$annotation_info)) {
                vars_to_select <- c(vars_to_select, ds_correlate$annotation_info$column)
            }
            vars_to_select <- c(vars_to_select, "symbol", "chr", "start", "end", "cor")
        }

        # get the indices into the annotype data
        annotations <- ds_correlate$annotations
        idxs <- match(names(pcor), annotations$annotation_id)

        correlations <-
            correlations %>%
            dplyr::inner_join(
                ds_correlate$annotations,
                by = c("annotation_id")
            ) %>%
            dplyr::select(vars_to_select)
    }

    if (valid(correlations)) {
        # don't return the NAs
        correlations <- correlations %>%
            dplyr::filter(!is.na(.data$cor))

        # attach some attributes
        attr(correlations, 'imputed_samples') <- samples_imputed
        attr(correlations, 'covar_formula') <- covar_formula
    }

    correlations
}

