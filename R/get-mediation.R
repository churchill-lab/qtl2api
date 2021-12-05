#' Perform mediation
#'
#' @param ds the dataset object
#' @param id the unique id in the dataset
#' @param marker_id marker identifier
#' @param ds_mediate the dataset object to mediate against
#'
#' @return A data.frame with the following columns depending on datatype:
#'         mRNA = gene_id, symbol, chr, pos, LOD
#'         protein = protein_id, gene_id, symbol, chr, pos, LOD
#'         phenotype = NONE
#'
#' @importFrom rlang .data
#' @export
get_mediation <- function(ds, id, marker_id, ds_mediate = NULL) {
    # get the data
    data <- get_data(ds)

    # get the dataset we are mediating against and the data
    ds_mediate <- nvl(ds_mediate, ds)
    data_mediate <- get_data(ds_mediate)

    # check if id exists
    idx <- which(colnames(data) == id)

    if (gtools::invalid(idx)) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # get the marker index and check it
    markers_cleaned <- markers %>% janitor::clean_names()
    mrkx <- which(markers_cleaned$marker_id == marker_id)

    if (gtools::invalid(mrkx)) {
        stop(sprintf("Cannot find marker '%s' in markers", marker_id))
    }

    # get the annotations
    if (tolower(ds_mediate$datatype) == "mrna") {
        # grab the annotations, create middle_point, and select what is needed
        annot <-
            ds_mediate$annot.mrna %>%
            janitor::clean_names() %>%
            dplyr::inner_join(
                tibble::enframe(colnames(data_mediate), name = NULL),
                by = c("gene_id" = "value")
            ) %>%
            dplyr::mutate(
                mid_point = (.data$start + .data$end) / 2
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
            ds_mediate$annot.protein %>%
            janitor::clean_names() %>%
            dplyr::inner_join(
                tibble::enframe(colnames(data_mediate), name = NULL),
                by = c("protein_id" = "value")
            ) %>%
            dplyr::mutate(
                mid_point = (.data$start + .data$end) / 2
            ) %>%
            dplyr::select(
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

    # get the covar data
    covar <- get_covar_matrix(ds_mediate)

    chrom <- as.character(markers[mrkx, "chr"])

    filtered_genoprobs <-
        genoprobs[[chrom]][rownames(data_mediate), , marker_id]

    intermediate::mediation.scan(
        target = data[, idx, drop = FALSE],
        mediator = data_mediate,
        annotation = annot,
        covar = covar,
        qtl.geno = filtered_genoprobs,
        verbose = FALSE
    )
}
