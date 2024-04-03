#' Create the covar.matrix element.
#'
#' @param dataset the dataset object (synchronized)
#' @param id the phenotype identifier, for phenotype datasets
#'
#' @return a named list list with 2 elements, covar_matrix and covar_formula
#'
#' @importFrom rlang .data
get_covar_matrix <- function(dataset, id = NULL) {
    # make sure samples and annotations are available
    ds <- synchronize_dataset(dataset)

    covar_formula <- NULL
    covar_matrix <- NULL

    if (is_phenotype(ds)) {
        # get the annot_phenotype row to get use_covar variable from the
        # annotations
        pheno <-
            ds$annot_phenotype %>%
            dplyr::filter(.data$data_name == id)

        if (qtl2api::invalid(pheno)) {
            stop(sprintf("Cannot find phenotype '%s' in dataset", id))
        }

        if (qtl2api::valid(pheno$use_covar)) {
            # create a string (model formula) from the use.covar column
            covar_formula <- paste0("~", gsub(":", "+", pheno$use_covar))
        }
    } else {
        if (qtl2api::valid(ds$covar_info)) {
            covar_formula <- paste0(ds$covar_info$sample_column, collapse="+")
            covar_formula <- paste0("~", covar_formula)
        }
    }

    if (qtl2api::valid(covar_formula)) {
        # get the sample id field
        sample_id_field <- ds$sample_id_field

        # convert samples to data.frame because QTL2 relies heavily
        # on rownames and colnames, rownames currently are or will
        # soon be deprecated in tibbles
        samples <- as.data.frame(ds$annot_samples)

        # create the model matrix, we use na.action = stats::na.pass so we can set
        # the rownames below.  We than use na.omit to filter down the data.
        covar_matrix <- stats::model.matrix.lm(
            stats::as.formula(covar_formula),
            data = samples,
            na.action = stats::na.pass
        )

        # drop the Intercept column
        covar_matrix <- covar_matrix[, -1, drop = FALSE]

        # drop the covar column if it has all identical values
        covar_matrix <- covar_matrix %>%
            tibble::as_tibble() %>%
            dplyr::select_if(function(col) length(unique(col))>1)

        # convert to a matrix and set the rownames so scan1 will work
        covar_matrix <- as.matrix(covar_matrix)

        rownames(covar_matrix) <-
            (samples %>% dplyr::select(dplyr::matches(sample_id_field)))[[1]]

        # do not need NA values
        covar_matrix <- stats::na.omit(covar_matrix)
    }

    list(
        covar_formula  = covar_formula,
        covar_matrix   = covar_matrix
    )
}
