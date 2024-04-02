#' Get analysis of the environment for debugging purposes.
#'
#' @return A list with various information regarding the datasets.
#' @export
get_analysis <- function() {

    n <- 0
    for (g in genoprobs){
        n <- n + dim(g)[3]
    }

    info <- list(
        genoprobs_samples = dim(g)[1],
        genoprobs_strains = dim(g)[2],
        genoprobs_markers = n,
        marker_num        = dim(markers)[1],
        marker_pos        = ifelse(all(markers$pos < 1000), 'DECIMAL', 'INTEGER')
    )

    datasets <- get_dataset_info()

    info[['datasets']] <- datasets$datasets
    info[['ensembl_version']] <- datasets$ensembl_version

    return(info)
}

