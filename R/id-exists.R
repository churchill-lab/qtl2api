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
            all_ids <- ds$annot_mrna |>
                dplyr::filter(gene_id == id) |>
                dplyr::select(
                    id      = gene_id,
                    gene_id = gene_id
                )
        } else if(tolower(ds$datatype) == 'protein') {
            all_ids <- ds$annot_protein |>
                dplyr::filter(gene_id == id) |>
                dplyr::select(
                    id         = protein_id,
                    gene_id    = gene_id,
                    protein_id = protein_id
                )
        } else if(tolower(ds$datatype) == 'protein_uniprot') {
            all_ids <- ds$annot_protein_uniprot |>
                dplyr::filter(gene_id == id) |>
                dplyr::select(
                    id         = uniprot_id,
                    gene_id    = gene_id,
                    protein_id = protein_id,
                    uniprot_id = uniprot_id
                )
        } else if(tolower(ds$datatype) == 'phos') {
            all_ids <- ds$annot_phos |>
                dplyr::filter(gene_id == id) |>
                dplyr::select(
                    id         = phos_id,
                    gene_id    = gene_id,
                    protein_id = protein_id,
                    uniprot_id = uniprot_id,
                    phos_id    = phos_id
                )
        } else if(tolower(ds$datatype) == 'ptm') {
            all_ids <- ds$annot_ptm |>
                dplyr::filter(gene_id == id) |>
                dplyr::select(
                    id         = ptm_id,
                    gene_id    = gene_id,
                    protein_id = protein_id,
                    uniprot_id = uniprot_id,
                    peptide_id = peptide_id,
                    ptm_id     = ptm_id
                )
        } else if(tolower(ds$datatype) == 'peptide') {
            all_ids <- ds$annot_peptide |>
                dplyr::filter(gene_id == id) |>
                dplyr::select(
                    id         = peptide_id,
                    gene_id    = gene_id,
                    protein_id = protein_id,
                    uniprot_id = uniprot_id,
                    peptide_id = peptide_id
                )
        } else if(tolower(ds$datatype) == 'phenotype') {
            all_ids <- ds$annot_phenotype |>
                dplyr::filter(data_name == id) |>
                dplyr::select(
                    id        = data_name,
                    data_name = data_name
                )
        }

        if (valid(all_ids)) {
            ret[[d]] <- all_ids
        }


        #if (id %in% all_ids) {
        #    ids <- ds$annot.ptm |> dplyr::filter(gene_id = 'ENSMUSG00000078515')



            #ret[[d]] <- list(
            #    dataset_id        = d,
            #    dataset_datatype  = ds$datatype,
            #    dataset_name      = ds$display_name,
            #    id                = id
            #)
        #}
    }

    if (length(ret) == 0) {
        ret <- NULL
    }

    ret
}

