#' Get the sample id field.
#'
#' @param dataset the dataset object
#' @return the name of the sample id field
get_sample_id_field <- function(dataset) {
    if (!gtools::invalid(dataset$is_synchronized)) {
        return(dataset$sample_id_field)
    }

    if (is_phenotype(dataset)) {
        annots_field <- grep("^annots?(\\.|_){1}pheno(type)?s?$",
                             names(dataset),
                             value = TRUE)

        sample_id_field <-
            dataset[[annots_field]] %>%
            janitor::clean_names() %>%
            dplyr::filter(.data$is_id == TRUE)

        if (NROW(sample_id_field) != 1) {
            stop("is_id == TRUE more than once in annot.phenotype")
        }

        sample_id_field <- sample_id_field$data_name
    } else {
        annots_field <- grep("^annots?(\\.|_){1}samples?$",
                             names(dataset),
                             value = TRUE)

        sample_id_field <- grep(
            "^mouse(\\.|_)?id$|^sample(\\.|_)?id$",
            colnames(dataset[[annots_field]]),
            value = TRUE,
            ignore.case = TRUE
        )

        if (length(sample_id_field) != 1) {
            stop("unable to find a sample id field in annot.samples")
        }
    }

    sample_id_field
}


#' Synchronize sample IDs between objects.
#'
#' If each object exists, we get the intersection of the sample IDs,
#' sort them and subset each object.  Sample IDs must be in rownames for
#' each parameter.
#'
#' @param pheno a matrix or data.frame containing the phenotypes
#' @param probs a numeric matrix containing the founder allele probabilities
#' @param expr a numeric matrix containing the gene expression data
#' @param covar a numeric matrix containing the mapping covariates
#'
#' @return list with four elements: `pheno`, `probs`, `expr` and `covar`. The
#' Sample IDs will all be in the same order.
#'
#' @export
synchronize_samples <- function(pheno, probs, expr, covar = NULL) {
    samples <- intersect(rownames(pheno), rownames(probs))
    samples <- intersect(samples, rownames(expr))

    if (!is.null(covar)) {
        samples <- intersect(samples, rownames(covar))
    }

    if (length(samples) == 0) {
        message("There are no samples in common")
    }

    samples <- sort(samples)
    pheno <- pheno[samples, , drop = FALSE]
    probs <- probs[samples, , drop = FALSE]
    expr <- expr[samples, , drop = FALSE]

    if (!missing(covar)) {
        covar <- covar[samples, , drop = FALSE]
    }

    list(
        pheno = pheno,
        probs = probs,
        expr = expr,
        covar = covar
    )
}
