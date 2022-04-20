#' Validate the dataset to make sure qtl2api can use it.
#'
#' @param dataset_a the dataset_id as a string identifier
#' @param dataset_b the dataset_id as a string identifier
#' @param extensive TRUE to perform scans on random elements in the data
#'
#' @export
compare_datasets <- function(dataset_a, dataset_b, extensive = FALSE) {

    # get the datasets
    if (is.character(dataset_a)) {
        dataset_a <- get_dataset_by_id(dataset_a)
    }
    if (is.character(dataset_b)) {
        dataset_b <- get_dataset_by_id(dataset_b)
    }

    # synchronize the datasets
    ds_a <- synchronize_dataset(dataset_a)
    ds_b <- synchronize_dataset(dataset_b)

    results <- list()

    if (ds_a$datatype == ds_b$datatype) {
        results$datatype <- 'EQUAL'
    } else {
        results$datatype <- 'NOT EQUAL'
    }

    if (ds_a$display_name == ds_b$display_name) {
        results$display_name <- 'EQUAL'
    } else {
        results$display_name <- 'NOT EQUAL'
    }

    results$annotations <- compare_annotations(ds_a, ds_b)
    results$covar_info <- compare_covar_info(ds_a, ds_b)
    results$samples <- compare_samples(ds_a, ds_b)

    results
}

compare_annotations <- function(ds_a, ds_b) {
    #
    # get annotation info for a
    #

     annots_a <- NULL
     ids_a <- NULL
     cols_a <- NULL

     if (ds_a$datatype == 'mrna') {
         annots_a <- ds_a$annot_mrna
         ids_a <- annots_a$gene_id
     } else if (ds_a$datatype == 'protein') {
         annots_a <- ds_a$annot_protein
         ids_a <- annots_a$protein_id
     } else if (ds_a$datatype == 'phos') {
         annots_a <- ds_a$annot_phos
         ids_a <- annots_a$phos_id
     } else if (ds_a$datatype == 'phenotype') {
         annots_a <- ds_a$annot_phenotype
         ids_a <- annots_a$data_name
     } else {
         warning('Unknown ds_a$datatype')
     }

     cols_a <- colnames(annots_a)

    #
    # get annotation info for b
    #

     annots_b <- NULL
     ids_b <- NULL
     cols_b <- NULL

     if (ds_b$datatype == 'mrna') {
         annots_b <- ds_b$annot_mrna
         ids_b <- annots_b$gene_id
     } else if (ds_b$datatype == 'protein') {
         annots_b <- ds_b$annot_protein
         ids_b <- annots_b$protein_id
     } else if (ds_b$datatype == 'phos') {
         annots_b <- ds_b$annot_phos
         ids_b <- annots_b$phos_id
     } else if (ds_b$datatype == 'phenotype') {
         annots_b <- ds_b$annot_phenotype
         ids_b <- annots_b$data_name
     } else {
         warning('Unknown ds_b$datatype')
     }

     cols_b <- colnames(annots_b)

     results <- list()

    #
    # get column (annotation) information
    #

    results$num_annots_a <- length(cols_a)
    results$num_annots_b <- length(cols_b)

    if (length(cols_a) == length(cols_a)) {
        results$num_annots <- 'EQUAL'
    } else {
        results$num_annots <- 'NOT EQUAL'
    }

    results$annots_only_a <- setdiff(cols_a, cols_b)
    results$annots_only_b <- setdiff(cols_b, cols_a)
    results$annots_both <- intersect(cols_a, cols_b)

    #
    # get row (id) information
    #

    results$num_ids_a <- length(ids_a)
    results$num_ids_b <- length(ids_b)

    if (length(ids_a) == length(ids_b)) {
        results$num_ids <- 'EQUAL'
    } else {
        results$num_ids <- 'NOT EQUAL'
    }

    results$ids_only_a <- setdiff(ids_a, ids_b)
    results$ids_only_b <- setdiff(ids_b, ids_a)
    results$ids_both <- intersect(ids_a, ids_b)

    #
    # return results
    #

    results
}


compare_covar_info <- function(ds_a, ds_b) {
     results <- list()

    #
    # get covar information for a
    #

    covar_info_a <- ds_a$covar_info
    covars_a <- NULL

    if (valid(covar_info_a)) {
        covars_a <- covar_info_a$sample_column
        results$num_covars_a <- length(covars_a)
    } else {
        results$num_covars_a <- 0
    }

    #
    # get covar information for b
    #

    covar_info_b <- ds_b$covar_info
    covars_b <- NULL

    if (valid(covar_info_b)) {
        covars_b <- covar_info_b$sample_column
        results$num_covars_b <- length(covars_b)
    } else {
        results$num_covars_b <- 0
    }

    #
    # get row (id) information
    #

    if (length(results$num_covars_a) == length(results$num_covars_b)) {
        results$num_covars <- 'EQUAL'
    } else {
        results$num_covars <- 'NOT EQUAL'
    }

    results$covars_only_a <- setdiff(covars_a, covars_b)
    results$covars_only_b <- setdiff(covars_b, covars_a)
    results$covars_both <- intersect(covars_a, covars_b)

    #
    # loop through the covar information
    #

    results$covar_matches <- list()

    if (valid(results$covars_both)) {
        for(i in 1:length(results$covars_both)) {
            sample_column <- results$covars_both[i]

            # rows from a
            sample_vals_a <-
                ds_a$annot_samples %>%
                dplyr::select(
                    id  = .data[[ds_a$sample_id_field]],
                    val = .data[[sample_column]]
                )

            # rows from b
            sample_vals_b <-
                ds_b$annot_samples %>%
                dplyr::select(
                    id  = .data[[ds_b$sample_id_field]],
                    val = .data[[sample_column]]
                )

            # combine
            matching_vals <-
                dplyr::inner_join(
                    sample_vals_a,
                    sample_vals_b,
                    by = c('id', 'val')
                )

            results$covar_matches[[sample_column]] <- list(
                samples_only_a = setdiff(sample_vals_a$id, matching_vals$id),
                samples_only_b = setdiff(sample_vals_b$id, matching_vals$id),
                samples_both   = matching_vals$id,
                num_samples    = length(matching_vals$id)
            )
        }
    }



    #
    # return results
    #

    results
}


compare_samples <- function(ds_a, ds_b) {
    #
    # get samples for a
    #

     samples_a <- ds_a$annot_samples
     ids_a <- samples_a[[ds_a$sample_id_field]]
     cols_a <- colnames(samples_a)

    #
    # get samples for b
    #

     samples_b <- ds_b$annot_samples
     ids_b <- samples_b[[ds_b$sample_id_field]]
     cols_b <- colnames(samples_b)

     results <- list()

    #
    # get column (annotation) information
    #

    results$num_annots_a <- length(cols_a)
    results$num_annots_b <- length(cols_b)

    if (length(cols_a) == length(cols_a)) {
        results$num_annots <- 'EQUAL'
    } else {
        results$num_annots_a <- 'NOT EQUAL'
    }

    results$annots_only_a <- setdiff(cols_a, cols_b)
    results$annots_only_b <- setdiff(cols_b, cols_a)
    results$annots_both <- intersect(cols_a, cols_b)

    #
    # get row (id) information
    #

    results$num_samples_a <- length(ids_a)
    results$num_samples_b <- length(ids_b)

    if (length(ids_a) == length(ids_b)) {
        results$num_samples <- 'EQUAL'
    } else {
        results$num_samples <- 'NOT EQUAL'
    }

    results$samples_only_a <- setdiff(ids_a, ids_b)
    results$samples_only_b <- setdiff(ids_b, ids_a)
    results$samples_both <- intersect(ids_a, ids_b)

    #
    # return results
    #

    results
}
