#' Check if id exists and has data in dataset.
#'
#' @param id the id to check
#'
#' @return `TRUE` if id contains data, `FALSE` otherwise
#' @export
id_exists <- function(id) {
    ret <- list()

    datasets <- utils::apropos('^dataset\\.*', ignore.case = TRUE)

    for (d in datasets) {
        ds <- synchronize_dataset(get_dataset_by_id(d))
        all_ids <- NULL

        if (tolower(ds$datatype) == 'mrna') {
            all_ids <- ds$annot_mrna$gene_id
        } else if(tolower(ds$datatype) == 'protein') {
            all_ids <- ds$annot_protein$gene_id
        } else if(tolower(ds$datatype) == 'protein_uniprot') {
            all_ids <- ds$annot_protein_uniprot$gene_id
        } else if(tolower(ds$datatype) == 'phos') {
            all_ids <- ds$annot_phos$gene_id
        } else if(tolower(ds$datatype) == 'ptm') {
            all_ids <- ds$annot_ptm$gene_id
        } else if(tolower(ds$datatype) == 'peptide') {
            all_ids <- ds$annot_peptide$gene_id
        } else if(tolower(ds$datatype) == 'phenotype') {
            all_ids <- ds$annot_phenotype$data_name
        }

        if (id %in% all_ids) {
            ret[[d]] <- list(
                dataset_id        = d,
                dataset_datatype  = ds$datatype,
                dataset_name      = ds$display_name,
                id                = id
            )
        }
    }

    if (length(ret) == 0) {
        ret <- NULL
    }

    ret
}

