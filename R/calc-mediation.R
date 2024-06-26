#' Perform mediation.
#'
#' @param dataset The dataset object.
#' @param id The unique id in the dataset.
#' @param marker_id The unique marker identifier.
#' @param dataset_mediate The dataset object to mediate against (defaults to
#'   dataset)
#'
#' @return A `tibble` with the following columns depending on datatype:
#' \itemize{
#'   \item mRNA - gene_id, symbol, chr, pos, LOD
#'   \item protein = protein_id, gene_id, symbol, chr, pos, LOD
#'   \item phenotype = NONE
#' }
#'
#' @importFrom rlang .data
#' @export
calc_mediation <- function(dataset, id, marker_id, dataset_mediate = NULL) {
    # make sure annotations, data, and samples are synchronized
    ds <- synchronize_dataset(dataset)

    # get the dataset we are mediating against
    if (is.null(dataset_mediate)) {
        ds_mediate <- ds
    } else {
        ds_mediate <- synchronize_dataset(dataset_mediate)
    }

    # check if id exists
    if (!any(id == colnames(ds$data))) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # get the marker index and check it
    markers_cleaned <-
        markers %>%
        dplyr::filter(!is.na(.data$pos)) %>%
        janitor::clean_names()

    mrkx <- which(markers_cleaned$marker_id == marker_id)

    if (!any(mrkx)) {
        stop(sprintf("Cannot find marker '%s' in markers", marker_id))
    }

    # get the annotations
    if (tolower(ds_mediate$datatype) == "mrna") {
        # grab the annotations, create middle_point, and select what is needed
        annot <-
            ds_mediate$annot_mrna %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names() %>%
            dplyr::inner_join(
                tibble::enframe(colnames(ds_mediate$data), name = NULL),
                by = c("gene_id" = "value")
            ) %>%
            dplyr::mutate(
                mid_point = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::select(
                gene_id      = .data$gene_id,
                symbol       = .data$symbol,
                chr          = .data$chr,
                middle_point = .data$mid_point
            )
    } else if (tolower(ds_mediate$datatype) == "protein") {
        # grab the annotations, create middle_point, and select what is needed
        annot <-
            ds_mediate$annot_protein %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names() %>%
            dplyr::inner_join(
                tibble::enframe(colnames(ds_mediate$data), name = NULL),
                by = c("protein_id" = "value")
            ) %>%
            dplyr::mutate(
                mid_point = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::select(
                protein_id   = .data$protein_id,
                gene_id      = .data$gene_id,
                symbol       = .data$symbol,
                chr          = .data$chr,
                middle_point = .data$mid_point
            )
    } else if (tolower(ds_mediate$datatype) == "protein_uniprot") {
        # grab the annotations, create middle_point, and select what is needed
        annot <-
            ds_mediate$annot_protein_uniprot %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names() %>%
            dplyr::inner_join(
                tibble::enframe(colnames(ds_mediate$data), name = NULL),
                by = c("uniprot_id" = "value")
            ) %>%
            dplyr::mutate(
                mid_point = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::select(
                uniprot_id   = .data$uniprot_id,
                protein_id   = .data$protein_id,
                gene_id      = .data$gene_id,
                symbol       = .data$symbol,
                chr          = .data$chr,
                middle_point = .data$mid_point
            )
    } else if (tolower(ds_mediate$datatype) == "phos") {
        # grab the annotations, create middle_point, and select what is needed
        annot <-
            ds_mediate$annot_phos %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names() %>%
            dplyr::inner_join(
                tibble::enframe(colnames(ds_mediate$data), name = NULL),
                by = c("phos_id" = "value")
            ) %>%
            dplyr::mutate(
                mid_point = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::select(
                phos_id      = .data$phos_id,
                protein_id   = .data$protein_id,
                gene_id      = .data$gene_id,
                symbol       = .data$symbol,
                chr          = .data$chr,
                middle_point = .data$mid_point
            )
    } else if (tolower(ds_mediate$datatype) == "ptm") {
        # grab the annotations, create middle_point, and select what is needed
        annot <-
            ds_mediate$annot_ptm %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names() %>%
            dplyr::inner_join(
                tibble::enframe(colnames(ds_mediate$data), name = NULL),
                by = c("ptm_id" = "value")
            ) %>%
            dplyr::mutate(
                mid_point = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::select(
                uniprot_id   = .data$uniprot_id,
                ptm_id       = .data$ptm_id,
                peptide_id   = .data$peptide_id,
                protein_id   = .data$protein_id,
                gene_id      = .data$gene_id,
                symbol       = .data$symbol,
                chr          = .data$chr,
                middle_point = .data$mid_point
            )
    } else if (tolower(ds_mediate$datatype) == "peptide") {
        # grab the annotations, create middle_point, and select what is needed
        annot <-
            ds_mediate$annot_peptide %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names() %>%
            dplyr::inner_join(
                tibble::enframe(colnames(ds_mediate$data), name = NULL),
                by = c("peptide_id" = "value")
            ) %>%
            dplyr::mutate(
                mid_point = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::select(
                uniprot_id   = .data$uniprot_id,
                peptide_id   = .data$peptide_id,
                protein_id   = .data$protein_id,
                gene_id      = .data$gene_id,
                symbol       = .data$symbol,
                chr          = .data$chr,
                middle_point = .data$mid_point
            )
    } else if (is_phenotype(ds_mediate)) {
        stop("dataset is a phenotype dataset and is not supported")
    } else {
        stop(sprintf("invalid dataset datatype: '%s'", ds_mediate$datatype))
    }

    # get the covar information
    covar_information <- get_covar_matrix(ds_mediate)

    chrom <- as.character(markers[mrkx, "chr"])

    samples <-
        ds_mediate$annot_samples %>%
        dplyr::select(dplyr::matches(ds_mediate$sample_id_field))

    # subset to the intersecting data
    sample_names <- intersect(samples[[1]], rownames(K[[chrom]]))
    sample_names <- intersect(sample_names, rownames(ds$data))
    sample_names <- intersect(sample_names, rownames(ds_mediate$data))

    filtered_genoprobs <-
        genoprobs[[chrom]][sample_names, , marker_id]

    data <- mediation.scan(
        target     = ds$data[sample_names, id, drop = FALSE],
        mediator   = ds_mediate$data[sample_names,],
        annotation = annot,
        covar      = covar_information$covar_matrix[sample_names, ],
        qtl.geno   = filtered_genoprobs,
        verbose    = FALSE
    )

    data$middle_point <- as.integer(data$middle_point)

    if (all(data$middle_point < 1000)) {
        data$middle_point <- as.integer(data$middle_point * 1000000)
    }

    attr(data, 'covar_formula') <- covar_information$covar_formula

    data
}

#' Mediation Scan.
#'
#' For a given QTL haplotype probabilities \code{qtl.geno} and target \code{target},
#' the function sequentially tries to add each column \code{m} of \code{mediator} matrix as a covariate
#' and calculates LOD statistic. The low LOD value indicates \code{qtl.geno} and
#' \code{target} are conditionally independent given \code{m},
#' i.e. \code{m} is a mediator of causal relationship from \code{qtl.geno} to \code{target}.
#'
#'
#' @param target A numeric vector with gene/protein expression
#' @param mediator A matrix, each column is one gene/protein's expression
#' @param annotation A data frame with mediators' annotation, must include columns "chr" and "pos"
#' @param qtl.geno A matrix, haplotype probabilities at QTL we try to mediate
#' @param covar A matrix with additive covariates
#' @param method A method to handle missing cases
#' @param verbose If TRUE display information about the progress
mediation.scan <- function(target,
                           mediator,
                           annotation,
                           qtl.geno,
                           covar=NULL,
                           method=c("double-lod-diff", "ignore", "lod-diff",
                                    "lod-ratio"),
                           verbose=TRUE) {

    # calculates log10-Likelihood of linear model y ~ 1 + X
    LL <- function(y, X) {
      -length(y)/2*log10(sum(qr.resid(qr(cbind(X,1)),y)^2))
    }

    # Synch sample IDs.
    tmp = synchronize_samples(
        pheno = target,
        probs = qtl.geno,
        expr = mediator,
        covar = covar
    )

    target   = tmp$pheno
    qtl.geno = tmp$probs
    mediator = tmp$expr
    covar    = tmp$covar
    rm(tmp)

    # check input
    stopifnot(NROW(target) == NROW(mediator))
    stopifnot(NROW(annotation) == NCOL(mediator))
    stopifnot(NROW(qtl.geno) == NROW(target))
    stopifnot(!any(is.na(qtl.geno)))
    stopifnot(all(is.numeric(target[,1])))
    stopifnot(all(is.numeric(mediator)))
    stopifnot(all(is.numeric(qtl.geno)))
    stopifnot(c("CHR", "MIDDLE_POINT") %in% toupper(names(annotation)))
    method = match.arg(method)

    if (!is.null(covar)) {
        stopifnot(NROW(target) == NROW(covar))
        stopifnot(!any(is.na(covar)))
        stopifnot(all(is.numeric(covar)))
    }

    # data preparation
    mediator <- cbind(mediator) # to ensure 'mediator' is a matrix
    N <- ncol(mediator) # number of points to scan
    if (is.null(covar)) covar <- cbind(rep(1, length(target))) # if no covariates, use just intercept
    LOD <- rep(NA, N) # prepare output

    if (method == "double-lod-diff") {
        no.na <- !is.na(target)
        LOD0 <-
            LL(target[no.na], cbind(covar, qtl.geno)[no.na,]) -
            LL(target[no.na], covar[no.na,])
    }

    # for-loop comparing M0: target~covar+mediator[,i] vs M1: target~covar+mediator[,i]+qtl.geno
    for (i in 1:N) {
        no.na <- !is.na(target) & !is.na(mediator[,i])
        loglik0 <- LL(target[no.na], cbind(covar[no.na,], mediator[no.na,i]))
        loglik1 <- LL(target[no.na], cbind(covar[no.na,], mediator[no.na,i], qtl.geno[no.na,]))

        if (method == "ignore" | (method == "double-lod-diff" & all(no.na))) {
            # "double-lod-diff" for no missing observation is identical to "ignore"
            LOD[i] <- loglik1 - loglik0
        } else {
            loglik2 <- LL(target[no.na,drop=FALSE], covar[no.na,,drop=FALSE])
            #loglik2 <- LL(target[no.na], covar[no.na,])
            loglik3 <- LL(target[no.na], cbind(covar[no.na,], qtl.geno[no.na,]))

            if (method == "lod-diff") {
                LOD[i] <- loglik3 - loglik2 - (loglik1-loglik0)
            } else if (method == "double-lod-diff") {
                LOD[i] <- LOD0 - (loglik3 - loglik2 - (loglik1-loglik0))
            } else if (method == "lod-ratio") {
                LOD[i] <- (10^loglik1-10^loglik0) / (10^loglik3 - 10^loglik2)
            }
        }
    }

    output <- annotation
    output$LOD <- LOD
    class(output) <- c("mediation", "data.frame")

    return(output)
}


