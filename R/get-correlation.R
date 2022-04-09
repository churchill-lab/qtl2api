#' Get the correlation.
#'
#' @param dataset the dataset object
#' @param id the unique id in the dataset
#' @param dataset_correlate the dataset to correlate against, dataset if INVALID
#' @param intcovar the interactive covariate
#' @param use_qr qr decomposition
#' @param backfill_NA if we need to impute, should we backfill the NAs
#'
#' @return a `tibble` with the correlation and annotations
#'
#' @export
get_correlation <- function(dataset, id, dataset_correlate = NULL,
                            intcovar = NULL, use_qr = TRUE,
                            backfill_NA = TRUE) {

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

    # make sure we have the same samples
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

    correlations <- NULL

    if (tolower(ds_correlate$datatype) == "mrna") {
        # get the indices into the annotype data
        annot_mrna <- ds_correlate$annot_mrna
        idxs <- match(names(pcor), annot_mrna$gene_id)

        correlations <- tibble::tibble(
            cor    = pcor,
            id     = names(pcor),
            symbol = annot_mrna$symbol[idxs],
            chr    = annot_mrna$chr[idxs],
            start  = as.integer(annot_mrna$start[idxs]),
            end    = as.integer(annot_mrna$end[idxs])
        )

        if (all(correlations$start < 1000)) {
            correlations$start <- as.integer(correlations$start * 1000000)
        }

        if (all(correlations$end < 1000)) {
            ret$end <- as.integer(correlations$end * 1000000)
        }
    } else if (tolower(ds_correlate$datatype) == "protein") {
        # get the indices into the annotype data
        annot_protein <- ds_correlate$annot_protein
        idxs <- match(names(pcor), annot_protein$protein_id)

        correlations <- tibble::tibble(
            cor     = pcor,
            id      = names(pcor),
            gene_id = annot_protein$gene_id[idxs],
            symbol  = annot_protein$symbol[idxs],
            chr     = annot_protein$chr[idxs],
            start   = as.integer(annot_protein$start[idxs]),
            end     = as.integer(annot_protein$end[idxs])
        )

        if (all(correlations$start < 1000)) {
            correlations$start <- as.integer(correlations$start * 1000000)
        }

        if (all(correlations$end < 1000)) {
            correlations$end <- as.integer(correlations$end * 1000000)
        }
    } else if (is_phenotype(ds_correlate)) {
        correlations <- tibble::tibble(
            cor = pcor,
            id = names(pcor)
        )
    }

    # don't return the NAs
    correlations <- correlations %>%
        dplyr::filter(!is.na(.data$cor))

    list(correlations    = correlations,
         imputed_samples = samples_imputed,
         covar_formula   = covar_formula)
}


#' Get the correlation data for plotting.
#'
#' @param dataset the dataset object
#' @param id the unique id in the dataset
#' @param dataset_correlate the dataset to correlate to
#' @param id_correlate the identifier from the correlate dataset
#' @param intcovar the interactive covariate
#' @param use_qr qr decomposition
#' @param backfill_NA if we need to impute, should we backfill the NAs
#'
#' @return a named `list` with the data to plot
#'
#' @export
get_correlation_plot_data <- function(dataset, id,
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
        # get the covar matrix
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
        dplyr::mutate(imputed = .data$sample_id %in% samples_imputed) %>%
        dplyr::filter(!is.na(.data$x)) %>%  # don't return the NAs
        dplyr::filter(!is.na(.data$y))      # don't return the NAs

    ret <- list(
        datatypes     = datatypes,
        data          = correlation_plot_data,
    )

    attr(ret, 'covar_formula') <- covar_information$covar_formula

    ret
}

#' Calculate the residual matrix
#'
#' @param variable_matrix The data  matrix for first set.
#' @param adjust_matrix The data matrix for the second set.
#' @param variables_interest List of variables of interest.
#' @param variables_compare List of variables to compare.
#' @param use_qr qr decomposition
#'
#' @return residual matrix
calc_residual_matrix <- function(variable_matrix,
                                 adjust_matrix,
                                 variables_interest,
                                 variables_compare,
                                 use_qr = TRUE) {

    # impute if necessary
    if(any(is.na(variable_matrix))) {
        variable_matrix <- missMDA::imputeFAMD(X = variable_matrix)$completeObs
    }

    # make sure we have the same samples
    samples <- intersect(rownames(variable_matrix), rownames(adjust_matrix))

    # combine the data
    data <- cbind(variable_matrix[samples, ], adjust_matrix[samples, ])
    data <- as.data.frame(data)

    # columns/variable names with ";" don't work well in formulas
    colnames(data) <-
        stringr::str_replace_all(colnames(data), ';', '_')
    variables_interest_fix <-
        stringr::str_replace_all(variables_interest, ';', '_')

    residual_matrix <- NULL

    if(use_qr) {
        formula_str <-
            paste("~ + 1 +", paste(variables_compare, collapse = " + "))

        X_0 <- stats::model.matrix.lm(stats::as.formula(formula_str),
                                      data = data,
                                      na.action = stats::na.exclude)

        X_0 <- X_0[samples, ]

        y_data <- data[samples, variables_interest_fix, drop = FALSE]
        colnames(y_data) <- variables_interest_fix

        qr_0 <- qr(X_0)

        residual_matrix <- sapply(seq(variables_interest_fix), function(i) {
            d <- y_data[, variables_interest_fix[i]]
            return(qr.resid(qr_0, d))
        }, simplify = TRUE)

        rownames(residual_matrix) <- samples
    }
    else{
        ## Way too slow, use QR trick
        residual_matrix <- sapply(seq(variables_interest_fix), function(i) {
            formula_str <- paste(
                variables_interest_fix[i],
                "~",
                paste(variables_compare, collapse = " + ")
            )

            fit <- stats::lm(
                stats::formula(formula_str),
                data = data,
                na.action = stats::na.exclude
            )

            return(fit$residuals[samples])
        }, simplify = TRUE)
    }

    colnames(residual_matrix) <- variables_interest

    residual_matrix
}
