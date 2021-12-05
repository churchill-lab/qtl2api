#' Get the correlation.
#'
#' @param ds the dataset object
#' @param id the unique id in the dataset
#' @param ds_correlate the dataset to correlate against, ds if INVALID
#' @param intcovar the interactive covariate
#' @param use_qr `TRUE` to use QR decomposition (FASTER)
#'
#' @return a `tibble` with the correlation and annotations
#'
#' @importFrom rlang .data
#' @export
get_correlation <- function(ds, id, ds_correlate = NULL, intcovar = NULL,
                            use_qr = TRUE) {
    # get the data
    data <- get_data(ds)

    # get the dataset we are correlating against and the data
    ds_correlate <- nvl(ds_correlate, ds)
    data_correlate <- get_data(ds_correlate)

    # check if id exists
    idx <- which(colnames(data) == id)

    if (gtools::invalid(idx)) {
        stop(sprintf("Cannot find data for '%s' in dataset", id))
    }

    # make sure we have the same samples
    samples <- intersect(rownames(data), rownames(data_correlate))

    if (length(samples) == 0) {
        stop("No matching samples for correlation")
    }

    data <- data[samples, ]
    data_correlate <- data_correlate[samples, ]

    # get the covar data
    covar <- get_covar_matrix(ds, id)

    if (gtools::invalid(intcovar)) {
        pcor <- stats::cor(data[, idx], data_correlate, use = "pair")
    } else {
        interactive_covariate <-
            colnames(covar)[grepl(intcovar, colnames(covar), ignore.case = T)]

        id_residual_matrix <-
            calc_residual_matrix(
                variable_matrix    = data,
                adjust_matrix      = covar,
                variables_interest = c(id),
                variables_compare  = interactive_covariate,
                use_qr             = use_qr
            )

        residual_matrix <-
            calc_residual_matrix(
                variable_matrix    = data_correlate,
                adjust_matrix      = covar,
                variables_interest = colnames(data_correlate),
                variables_compare  = interactive_covariate,
                use_qr             = use_qr
            )

        resid_samples <- intersect(
            rownames(id_residual_matrix),
            rownames(residual_matrix)
        )

        pcor <- stats::cor(
            id_residual_matrix[resid_samples, ],
            residual_matrix[resid_samples, ],
            use = "pair"
        )
    }

    # reorder in decreasing order
    pcor <- pcor[1, order(abs(pcor), decreasing = TRUE)]

    ret <- NULL

    if (tolower(ds_correlate$datatype) == "mrna") {
        # get the indices into the annotype data
        annot_mrna <- ds_correlate$annot.mrna %>% janitor::clean_names()
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
        annot_protein <- ds_correlate$annot.protein %>% janitor::clean_names()
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
#' @param ds the dataset object
#' @param id the unique id in the dataset
#' @param ds_correlate the dataset to correlate to
#' @param id_correlate the identifier from the correlate dataset
#' @param intcovar the interactive covariate
#'
#' @return a named `list` with the data to plot
#'
#' @importFrom rlang .data
#' @export
get_correlation_plot_data <- function(ds, id,
                                      ds_correlate, id_correlate,
                                      intcovar = NULL) {
    # get the data
    data <- get_data(ds)

    # get the correlation dataset and data
    ds_correlate <- nvl(ds_correlate, ds)
    data_correlate <- get_data(ds_correlate)

    # make sure we have the same samples
    samples <- intersect(rownames(data), rownames(data_correlate))

    if (length(samples) == 0) {
        stop("No matching samples to correlate")
    }

    data <- data[samples, ]
    data_correlate <- data_correlate[samples, ]

    # check if id exists and get the index
    idx <- which(colnames(data) == id)
    idx_correlate <- which(colnames(data_correlate) == id_correlate)

    if (gtools::invalid(idx)) {
        stop(sprintf("Cannot find data for '%s' in dataset", id))
    } else if (gtools::invalid(idx_correlate)) {
        stop(sprintf(
            "Cannot find data for '%s' in dataset (correlate)",
            id_correlate
        ))
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
                variables_compare  = interactive_covariate,
                use_qr             = TRUE
            )

        data_correlate <-
            calc_residual_matrix(
                variable_matrix    = data_correlate,
                adjust_matrix      = covar,
                variables_interest = colnames(data_correlate),
                variables_compare  = interactive_covariate,
                use_qr             = TRUE
            )

        x <- data[, 1]
        y <- data_correlate[, idx_correlate]
    } else {
        x <- data[, idx]
        y <- data_correlate[, idx_correlate]
    }

    # get the sample id field
    sample_id_field <- get_sample_id_field(ds)

    # get the intersecting samples and indices
    samples <- intersect(rownames(data), rownames(data_correlate))
    samples_idx <- which(ds$annot.samples[[sample_id_field]] %in% samples)

    # get the covar factors and their data levels
    sample_info <- list()
    datatypes <- list()
    covar_info <- ds$covar.info %>% janitor::clean_names()

    for (s in covar_info$sample_column) {
        stopifnot(!is.null(ds$annot.samples[[s]]))
        sample_info[[toString(s)]] <- ds$annot.samples[samples_idx, ][[s]]

        if (is.factor(ds$annot.samples[[s]])) {
            datatypes[[toString(s)]] <-
                gtools::mixedsort(levels(ds$annot.samples[[s]]))
        } else {
            datatypes[[toString(s)]] <-
                gtools::mixedsort(unique(ds$annot.samples[[s]]))
        }
    }

    correlation_plot_data <-
        tibble::as_tibble(data.frame(
            sample_id = rownames(data),
            x = x,
            y = y,
            sample_info,
            stringsAsFactors = FALSE
        ))

    # TODO: should we add id to to the datset object (fix_environemnt)
    list(
        #dataset           = nvl(ds$id, ds$display.name),
        id                = id,
        #dataset.correlate = nvl(ds_correlate$id, ds_correlate$display.name),
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
#' @param use_qr TRUE to use QR decomposition (FASTER).
#' @param impute TRUE to impute NA data.
#'
#' @return residual matrix
#'
calc_residual_matrix <- function(variable_matrix,
                                 adjust_matrix,
                                 variables_interest,
                                 variables_compare,
                                 use_qr = TRUE,
                                 impute = TRUE) {
    # make sure we have the same samples
    samples <- intersect(rownames(variable_matrix), rownames(adjust_matrix))

    # combine the data
    data <- cbind(variable_matrix[samples, ], adjust_matrix[samples, ])
    data <- as.data.frame(data)

    if (use_qr) {
        # Fast, but can't handle NAs in y??
        formula_str <-
            paste("~ + 1 +", paste(variables_compare, collapse = " + "))

        X_0 <- stats::model.matrix.lm(
            stats::as.formula(formula_str),
            data = data,
            na.action = stats::na.pass
        )

        X_0 <- X_0[samples, ]

        y_data <- data[samples, variables_interest, drop = FALSE]
        colnames(y_data) <- variables_interest

        if (impute) {
            if (any(is.na(X_0))) {
                X_0 <- missMDA::imputeFAMD(X = X_0)$completeObs
            }

            if (any(is.na(y_data))) {
                y_data <- missMDA::imputeFAMD(X = y_data)$completeObs
            }
        }

        qr_0 <- qr(X_0)

        residual_matrix <- sapply(seq(variables_interest), function(i) {
            d <- y_data[, variables_interest[i]]
            return(qr.resid(qr_0, d))
        }, simplify = TRUE)

        rownames(residual_matrix) <- samples
    } else {
        ## Way too slow, use QR trick
        residual_matrix <- sapply(seq(variables_interest), function(i) {
            formula_str <- paste(
                variables_interest[i],
                "~",
                paste(variables_compare, collapse = " + ")
            )
            fit <- stats::lm(stats::formula(formula_str), data = data)
            return(fit$residuals[samples])
        }, simplify = TRUE)
    }

    colnames(residual_matrix) <- variables_interest
    residual_matrix
}
