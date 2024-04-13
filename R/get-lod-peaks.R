#' Get the LOD peaks for a dataset.
#'
#' By default, this will get all LOD peaks stored in the `dataset$lod_peaks`
#' `list`.  If `intcovar` is specified, only the peaks for that are returned.
#'
#' @param ds The dataset object to get the LOD peaks for.
#'
#' @return A named `list` of all the LOD peaks. Each element in the `list` is a
#' `tibble` with the following columns depending upon `dataset$datatype`:
#' \itemize{
#'   \item mRNA - `marker,chr,pos,gene_id,symbol,gene_chrom,middle,lod`
#'   \item protein - `marker,chr,pos,protein_id,gene_id,symbol,gene_chrom,middle,lod`
#'   \item phenotype - `marker,chr,pos,data_name,short_name,description,lod`
#' }
#'
#' If the allele effects are stored, values `A-H` will also be included.
#' @export
get_lod_peaks <- function(ds) {
    # these functions should not be synchronized
    if (valid(ds$is_synchronized)) {
        stop("dataset should not be synchronized")
    }

    ds_synch <- synchronize_dataset(ds)

    peaks <- list()

    # get the additive LOD peaks
    peaks_additive <- get_lod_peaks_by_covar(ds)

    if (!is.null(peaks_additive)) {
        peaks <- list(additive = peaks_additive)
    }

    # get the rest
    for (i in seq(nrow(ds_synch$covar_info))) {
        cov_inf <- ds_synch$covar_info[i, ]

        if (cov_inf$interactive == TRUE) {
            peaks[[cov_inf$sample_column]] <-
                get_lod_peaks_by_covar(ds, cov_inf$sample_column)
        }
    }

    if (length(peaks) == 0) {
        peaks <- NULL
    }

    peaks
}
