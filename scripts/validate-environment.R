#' Validate the environment for defined objects.
#'
#' @param extensive TRUE to perform scans on random elements in the data
#'
#' @export
validate_environment <- function(by_api = FALSE) {
    cat("\n")
    cat("qtl2api tries to use gene_id/gene.id, protein_id/protein.id,\n")
    cat("data_name/data.name along as trying to figure out annotations for\n")
    cat("mouse_id/mouse.id/sample_id/sample.id\n")
    cat("\n")
    cat("Internally, qtl2api uses underscore '_' for variable names and\n")
    cat("column names\n")
    cat("\n")
    cat("If an error is detected, this method will use the '_' version\n")
    cat("\n")


    # grab the datasets in the environment
    datasets <- grep('^dataset*', utils::apropos('dataset\\.'), value=TRUE)

    if (invalid(datasets)) {
        message("ERROR   : No datasets found!")
        return()
    }

    ensembl_version <-
        utils::apropos("^ensembl(\\.|_){1}version$", ignore.case = TRUE)

    if (length(ensembl_version) == 0) {
        message("ERROR   : ensembl_version does not exist")
    }

    # Check genoprobs and K.
    if(length(genoprobs) != length(K)) {
        message("ERROR   : genoprobs (", length(genoprobs), ") and K (", length(K), ") do not have the same length.")
    } else {
        if(any(names(genoprobs) != names(K))) {
            message("ERROR   : names of genoprobs and K do not match.")
        }

        rownames.eq <- mapply(function(x, y) { all(rownames(x) == rownames(y)) }, genoprobs, K)
        if(any(rownames.eq == FALSE)) {
            message("ERROR   : sample IDs do not match between genoprobs and K.")
        }
    }

    # Check Marker IDs for genoprobs and map
    if(length(genoprobs) != length(map)) {
        message("ERROR   : genoprobs (", length(genoprobs), ") and map (", length(map), ") do not have the same length.")
    } else {
        rownames.eq <- mapply(function(x, y) { all(dimnames(x)[[3]] == names(y)) }, genoprobs, map)
        if(any(rownames.eq == FALSE)) {
            message("ERROR   : marker names do not match between genoprobs and map")
        }
    }

    # Check dimensions of markers and map.
    map.length = sum(sapply(map, length))
    if(map.length != nrow(markers)) {
        message("ERROR   : number of rows in markers (", nrow(markers), ") does not equal the number of markers in map (", map.length, ")")
    }

    for (ds in datasets) {
        cat("\nSTATUS  : Checking dataset: ", ds, " ...\n")
        validate_dataset(ds, by_api)
    }
}
