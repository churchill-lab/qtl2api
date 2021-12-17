#' Get the correlation.
#'
#' @param dataset the dataset object
#' @param id the unique id in the dataset
#' @param dataset_correlate the dataset to correlate against, dataset if INVALID
#' @param intcovar the interactive covariate
#'
#' @return a `tibble` with the correlation and annotations
#'
#' @importFrom rlang .data
#' @export
get_correlation <- function(dataset, id, dataset_correlate = NULL,
                            intcovar = NULL) {
    # make sure annotations, data, and samples are synchronized
    ds <- synchronize_dataset(dataset)

    # get the dataset we are correlating against
    if (gtools::invalid(dataset_correlate)) {
        ds_correlate <- ds
    } else {
        ds_correlate <- synchronize_dataset(dataset_correlate)
    }

    # check if id exists
    if (id %not in% colnames(ds$data)) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # make sure we have the same samples
    samples <- intersect(rownames(ds$data), rownames(ds_correlate$data))

    if (length(samples) == 0) {
        stop("No matching samples for correlation")
    }

    data <- ds$data[samples, ]
    data_correlate <- ds_correlate$data[samples, ]

    # get the covar data
    covar <- get_covar_matrix(ds, id)

    if (gtools::invalid(intcovar)) {
        pcor <- stats::cor(data[, id], data_correlate, use = "pair")
    } else {
        interactive_covariate <-
            colnames(covar)[grepl(intcovar, colnames(covar), ignore.case = T)]

        id_residual_matrix <-
            calc_residual_matrix(
                variable_matrix    = data,
                adjust_matrix      = covar,
                variables_interest = c(id),
                variables_compare  = interactive_covariate
            )

        residual_matrix <-
            calc_residual_matrix(
                variable_matrix    = data_correlate,
                adjust_matrix      = covar,
                variables_interest = colnames(data_correlate),
                variables_compare  = interactive_covariate
            )

        pcor <- stats::cor(
            id_residual_matrix,
            residual_matrix,
            use = "pair"
        )
    }

    # reorder in decreasing order
    pcor <- pcor[1, order(abs(pcor), decreasing = TRUE)]

    ret <- NULL

    if (tolower(ds_correlate$datatype) == "mrna") {
        # get the indices into the annotype data
        annot_mrna <- ds_correlate$annot_mrna
        idxs <- match(names(pcor), annot_mrna$gene_id)

        ret <- tibble::tibble(
            cor    = pcor,
            id     = names(pcor),
            symbol = annot_mrna$symbol[idxs],
            chr    = annot_mrna$chr[idxs],
            start  = annot_mrna$start[idxs],
            end    = annot_mrna$end[idxs]
        )
    } else if (tolower(ds_correlate$datatype) == "protein") {
        # get the indices into the annotype data
        annot_protein <- ds_correlate$annot_protein
        idxs <- match(names(pcor), annot_protein$protein_id)

        ret <- tibble::tibble(
            cor     = pcor,
            id      = names(pcor),
            gene_id = annot_protein$gene_id[idxs],
            symbol  = annot_protein$symbol[idxs],
            chr     = annot_protein$chr[idxs],
            start   = annot_protein$start[idxs],
            end     = annot_protein$end[idxs]
        )
    } else if (is_phenotype(ds_correlate)) {
        ret <- tibble::tibble(
            cor = pcor,
            id = names(pcor)
        )
    }

    ret
}


#' Get the correlation data for plotting.
#'
#' @param dataset the dataset object
#' @param id the unique id in the dataset
#' @param dataset_correlate the dataset to correlate to
#' @param id_correlate the identifier from the correlate dataset
#' @param intcovar the interactive covariate
#'
#' @return a named `list` with the data to plot
#'
#' @importFrom rlang .data
#' @export
get_correlation_plot_data <- function(dataset, id,
                                      dataset_correlate, id_correlate,
                                      intcovar = NULL) {
    # make sure samples and annotations are available
    ds <- synchronize_dataset(dataset)

    # get the dataset we are correlating against
    if (gtools::invalid(dataset_correlate)) {
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
    if (id %not in% colnames(ds$data)) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    } else if (id_correlate %not in% colnames(ds_correlate$data)) {
        stop(sprintf("Cannot find id '%s' in dataset", id_correlate))
    }

    # get the covar matrix
    covar <- get_covar_matrix(ds, id)

    if (!gtools::invalid(intcovar)) {
        interactive_covariate <-
            colnames(covar)[grepl(intcovar, colnames(covar), ignore.case = T)]

        data <-
            calc_residual_matrix(
                variable_matrix    = data,
                adjust_matrix      = covar,
                variables_interest = c(id),
                variables_compare  = interactive_covariate
            )

        data_correlate <-
            calc_residual_matrix(
                variable_matrix    = data_correlate,
                adjust_matrix      = covar,
                variables_interest = colnames(data_correlate),
                variables_compare  = interactive_covariate
            )

        x <- data[, 1]
        y <- data_correlate[, id_correlate]
    } else {
        x <- data[, id]
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
                sample_id = rownames(data),
                x         = x,
                y         = y,
                sample_info,
                stringsAsFactors = FALSE
        ))

    # TODO: should we add id to to the dataset object (fix_environemnt)
    list(
        #dataset           = nvl(ds$id, ds$display.name),
        #dataset.correlate = nvl(ds_correlate$id, ds_correlate$display.name),
        id                = id,
        id_correlate      = id_correlate,
        datatypes         = datatypes,
        data              = correlation_plot_data
    )
}

#' Calculate the residual matrix
#'
#' @param variable_matrix The data  matrix for first set.
#' @param adjust_matrix The data matrix for the second set.
#' @param variables_interest List of variables of interest.
#' @param variables_compare List of variables to compare.
#'
#' @return residual matrix
#'
calc_residual_matrix <- function(variable_matrix,
                                 adjust_matrix,
                                 variables_interest,
                                 variables_compare) {
    # make sure we have the same samples
    samples <- intersect(rownames(variable_matrix), rownames(adjust_matrix))

    # combine the data
    data <- cbind(variable_matrix[samples, ], adjust_matrix[samples, ])
    data <- as.data.frame(data)

    residual_matrix <- NULL

    tryCatch(
        {
            # Fast, but can't handle NAs in y??
            formula_str <-
                paste("~ + 1 +", paste(variables_compare, collapse = " + "))

            X_0 <- stats::model.matrix.lm(
                stats::as.formula(formula_str),
                data = data,
                na.action = stats::na.exclude
            )

            X_0 <- X_0[samples, ]

            y_data <- data[samples, variables_interest, drop = FALSE]

            colnames(y_data) <- variables_interest

            qr_0 <- qr(X_0)

            residual_matrix <- sapply(seq(variables_interest), function(i) {
                d <- y_data[, variables_interest[i]]
                return(qr.resid(qr_0, d))
            }, simplify = TRUE)

            rownames(residual_matrix) <- samples
            colnames(residual_matrix) <- variables_interest
        },
        error = function(cond) {
        },
        warning = function(cond) {
        },
        finally = {
        }
    )

    if (is.null(residual_matrix)) {
        ## Way too slow
        residual_matrix <- sapply(seq(variables_interest), function(i) {
            formula_str <- paste(
                variables_interest[i],
                "~",
                paste(variables_compare, collapse = " + ")
            )
            fit <- stats::lm(stats::formula(formula_str), data = data, na.action = na.exclude)
            return(fit$residuals[samples])
        }, simplify = TRUE)
        colnames(residual_matrix) <- variables_interest
    }

    residual_matrix
}
