#' Validate the lod.peaks to make sure qtl2api can use it.
#'
#' @param dataset the dataset as a string identifier
#'
#' @export
validate_lod_peaks <- function(dataset) {
    cat("STATUS  : Checking lod_peaks\n")

    # get the dataset
    if (is.character(dataset)) {
        ds_orig <- get_dataset_by_id(dataset)
    } else {
        ds_orig <- dataset
    }

    ds <- synchronize_dataset(ds_orig)

    annots_field <- grep("^lod(\\.|_){1}peaks$",
                         names(ds_orig),
                         value = TRUE)

    if (length(annots_field) == 0) {
        message("WARNING : lod_peaks not found")
        return()
    }

    if (!is.list(ds_orig[[annots_field]])) {
        message("ERROR   : lod_peaks should be a list, but found ", class(ds_orig[[annots_field]]))
    }

    if ('additive' %not in% names(ds_orig[[annots_field]])) {
        message("ERROR   : additive should be an element in lod_peaks")
    }

    # get the rest
    for (i in seq(nrow(ds$covar_info))) {
        cov_inf <- ds$covar_info[i, ]

        # only look at interactive peaks
        if (cov_inf$interactive) {
            annots_field_peaks <- grep("^lod(\\.|_){1}peaks?$",
                                       names(ds_orig),
                                       value = TRUE)

            peaks <- ds_orig[[annots_field_peaks]][[cov_inf$lod_peaks]]

            if (invalid(peaks)) {
                message("ERROR   : Unable to find lod_peaks '", cov_inf$lod_peaks, "'")
            } else {

                peaks <- peaks %>%
                    janitor::clean_names()

                if (tolower(ds$datatype) == "mrna") {
                    if (length(setdiff(peaks$gene_id, ds$annot_mrna$gene_id))) {
                        cat("WARNING : not all lod_peaks$gene_id are in annot_mrna$gene_id\n")
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        cat("WARNING : not all lod_peaks$marker_id are in markers\n")
                    }
                } else if (tolower(ds$datatype) == "protein") {
                    if (length(setdiff(peaks$protein_id, ds$annot_protein$protein_id))) {
                        cat("WARNING : not all lod_peaks$protein_id are in annot_protein$protein_id\n")
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        cat("WARNING : not all lod_peaks$marker_id are in markers\n")
                    }
                } else if (tolower(ds$datatype) == "protein_uniprot") {
                    if (length(setdiff(peaks$uniprot_id, ds$annot_protein_uniprot$uniprot_id))) {
                        cat("WARNING : not all lod_peaks$uniprot_id are in annot_protein_uniprot$uniprot_id\n")
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        cat("WARNING : not all lod_peaks$marker_id are in markers\n")
                    }
                } else if (tolower(ds$datatype) == "phos") {
                    if (length(setdiff(peaks$phos_id, ds$annot_phos$phos_id))) {
                        cat("WARNING : not all lod_peaks$phos_id are in annot_phos$phos_id\n")
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        cat("WARNING : not all lod_peaks$marker_id are in markers\n")
                    }
                } else if (tolower(ds$datatype) == "ptm") {
                    if (length(setdiff(peaks$ptm_id, ds$annot_ptm$ptm_id))) {
                        cat("WARNING : not all lod_peaks$ptm_id are in annot_ptm$ptm_id\n")
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        cat("WARNING : not all lod_peaks$marker_id are in markers\n")
                    }
                } else if (tolower(ds$datatype) == "peptide") {
                    if (length(setdiff(peaks$peptide_id, ds$annot_peptide$peptide_id))) {
                        cat("WARNING : not all lod_peaks$peptide_id are in annot_peptide$peptide_id\n")
                    }

                    if (length(setdiff(peaks$marker_id, markers$marker.id))) {
                        cat("WARNING : not all lod_peaks$marker_id are in markers\n")
                    }
                } else {
                    if (length(setdiff(peaks$data_name, ds$annot_phenotype$data_name))) {
                        cat("WARNING : not all lod_peaks['", cov_inf$sample_column, "']$data_name are in annot_phenotype$data_name\n")
                    }
                }
            }
        }
    }
}
