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
#'
#' @return A `tibble` with the correlation and annotations.
#'
#' @export
calc_correlation <- function(dataset, id, dataset_correlate = NULL,
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

    # attach some annotations to the correlation data

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
    } else if (tolower(ds_correlate$datatype) == "protein_uniprot") {
        # get the indices into the annotype data
        annot_protein_uniprot <- ds_correlate$annot_protein_uniprot
        idxs <- match(names(pcor), annot_protein_uniprot$uniprot_id)

        correlations <- tibble::tibble(
            cor        = pcor,
            id         = names(pcor),
            protein_id = annot_protein_uniprot$protein_id[idxs],
            gene_id    = annot_protein_uniprot$gene_id[idxs],
            symbol     = annot_protein_uniprot$symbol[idxs],
            chr        = annot_protein_uniprot$chr[idxs],
            start      = as.integer(annot_protein_uniprot$start[idxs]),
            end        = as.integer(annot_protein_uniprot$end[idxs])
        )
    } else if (tolower(ds_correlate$datatype) == "phos") {
        # get the indices into the annotype data
        annot_phos <- ds_correlate$annot_phos
        idxs <- match(names(pcor), annot_phos$phos_id)

        correlations <- tibble::tibble(
            cor        = pcor,
            id         = names(pcor),
            protein_id = annot_phos$protein_id[idxs],
            gene_id    = annot_phos$gene_id[idxs],
            symbol     = annot_phos$symbol[idxs],
            chr        = annot_phos$chr[idxs],
            start      = as.integer(annot_phos$start[idxs]),
            end        = as.integer(annot_phos$end[idxs])
        )
    } else if (tolower(ds_correlate$datatype) == "ptm") {
        # get the indices into the annotype data
        annot_ptm <- ds_correlate$annot_ptm
        idxs <- match(names(pcor), annot_ptm$ptm_id)

        correlations <- tibble::tibble(
            cor        = pcor,
            id         = names(pcor),
            ptm_id     = annot_ptm$ptm_id[idxs],
            peptide_id = annot_ptm$peptide_id[idxs],
            protein_id = annot_ptm$protein_id[idxs],
            gene_id    = annot_ptm$gene_id[idxs],
            symbol     = annot_ptm$symbol[idxs],
            uniprot_id = annot_ptm$uniprot_id[idxs],
            chr        = annot_ptm$chr[idxs],
            start      = as.integer(nvl(annot_ptm$start[idxs], 0)),
            end        = as.integer(nvl(annot_ptm$end[idxs], 0))
        )
    } else if (tolower(ds_correlate$datatype) == "peptide") {
        # get the indices into the annotype data
        annot_peptide <- ds_correlate$annot_peptide
        idxs <- match(names(pcor), annot_peptide$ptm_id)

        correlations <- tibble::tibble(
            cor        = pcor,
            id         = names(pcor),
            peptide_id = annot_peptide$peptide_id[idxs],
            protein_id = annot_peptide$protein_id[idxs],
            gene_id    = annot_peptide$gene_id[idxs],
            symbol     = annot_peptide$symbol[idxs],
            uniprot_id = annot_peptide$uniprot_id[idxs],
            chr        = annot_peptide$chr[idxs],
            start      = as.integer(nvl(annot_peptide$start[idxs], 0)),
            end        = as.integer(nvl(annot_peptide$end[idxs], 0))
        )
    } else if (is_phenotype(ds_correlate)) {
        correlations <- tibble::tibble(
            cor = pcor,
            id = names(pcor)
        )
    }

    if (!is_phenotype(ds_correlate)) {
        if (all(correlations$start < 1000)) {
            correlations$start <- as.integer(correlations$start * 1000000)
        }

        if (all(correlations$end < 1000)) {
            correlations$end <- as.integer(correlations$end * 1000000)
        }
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

