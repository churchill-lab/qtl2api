#' Calculate the residual matrix.
#'
#' @param variable_matrix The data `matrix` for first set.
#' @param adjust_matrix The data `matrix` for the second set.
#' @param variables_interest `list` of variables of interest.
#' @param variables_compare `list` of variables to compare.
#' @param use_qr `TRUE` to use qr decomposition.
#'
#' @return The residual `matrix`.
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
    data <- cbind(
        variable_matrix[samples, , drop = FALSE],
        adjust_matrix[samples, , drop = FALSE]
    )
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
