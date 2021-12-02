
#' Get the SNP association mapping
#'
#' @param ds the dataset object
#' @param id the unique id in the dataset
#' @param chrom The chromosome.
#' @param intcovar the interactive covariate
#' @param chrom The chromosome.
#' @param location location on chromosome in base pairs
#' @param db_file full path to the sqlite database file
#' @param window_size the size of the window to scan before and after location
#' @param cores number of cores to use (0=ALL)
#'
#' @return a `data.frame` with the following columns: snp, chr, pos, alleles,
#'   sdp, ensembl_gene, csq, index, interval, on_map, lod.
#'
#' @importFrom rlang .data
#' @export
get_snp_assoc_mapping <- function(ds, id, intcovar, chrom, location, db_file,
                                  window_size = 500000, cores = 0) {
    # get the data
    data <- get_data(ds)

    # check if id exists
    idx <- which(colnames(data) == id)

    if (gtools::invalid(idx)) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    # make sure location is valid
    if (nvl_int(location, -1) == -1) {
        stop(paste0("Specified location is invalid: ", location))
    }

    position <- nvl_int(location, -1) / 1000000.0

    if (nvl_int(window_size, -1) == -1) {
        stop(paste0("Specified window_size is invalid: ", window_size))
    }

    window_size <- nvl_int(window_size, -1)
    window_length <- window_size / 1000000.0
    window_range_start <- max(position - window_length, 0)
    window_range_end <- position + window_length

    # make sure nCores is appropriate
    num_cores <- nvl_int(cores, 0)

    have_snps <- FALSE
    tries <- 1

    # extract SNPs from the database, we allow up to 10 tries to
    # find SNPs in a window
    while (!have_snps) {
        db_snps <- DBI::dbConnect(RSQLite::SQLite(), db_file)
        #db_snps <- dplyr::tbl(db_file, create = FALSE)
        window_snps <-
            dplyr::tbl(db_snps, "snps") %>%
            dplyr::filter(
                .data$chr     == chrom,
                .data$pos_Mbp >= window_range_start,
                .data$pos_Mbp <= window_range_end
            ) %>%
            dplyr::arrange(.data$pos_Mbp) %>%
            dplyr::collect(n = Inf)

        if (NROW(window_snps) > 0) {
            have_snps <- TRUE
        } else {
            window_size <- window_size + nvl_int(window_size, -1)
            window_length <- window_size / 1000000.0
            window_range_start <- max(position - window_length, 0)
            window_range_end <- position + window_length
            tries <- tries + 1

            if (tries > 10) {
                stop(sprintf(
                    "Cannot find any snps in region: %s:%f-%f",
                    chrom, window_range_start, window_range_end
                ))
            }
        }
    }

    colnames(window_snps)[c(1, 3)] <- c("snp", "pos")
    window_snps <- qtl2::index_snps(map = map, window_snps)

    # get the covar data
    covar <- get_covar_matrix(ds, id)

    # convert allele probs to SNP probs
    snp_prob <- qtl2::genoprob_to_snpprob(genoprobs, window_snps)

    # perform the scan using QTL2,
    # - addcovar should always be ALL covars
    # - intcovar should be just the interactive covariate column
    out_snps <- qtl2::scan1(
        pheno     = data[, idx, drop = FALSE],
        kinship   = K[[chrom]],
        genoprobs = snp_prob,
        addcovar  = covar,
        cores     = num_cores
    )

    map_tmp <- snpinfo_to_map(window_snps)
    tmp <- expand_snp_results(out_snps, map_tmp, window_snps)

    ret <- window_snps
    ret$lod <- tmp$lod[, 1]

    # set the interactive_covariates, to be used in scan1
    # as scan1(intcovar=interactive.covariate)
    interactive_covariate <- NULL

    if (!gtools::invalid(intcovar)) {
        if (intcovar %not in% ds$covar.info$sample_column) {
            stop(sprintf(
                "intcovar '%s' not found in covar.info in dataset",
                intcovar
            ))
        }

        # grabbing all the columns from covar (covar.matrix) that
        # match, i.e., "batch" will match "batch2", "BATCH3", etc
        interactive_covariate <-
            covar[, which(grepl(intcovar, colnames(covar), ignore.case = T))]

        # perform the scan using QTL2,
        # - addcovar should always be ALL covars
        # - intcovar should be just the interactive covariate column
        out_snps <- qtl2::scan1(
            pheno     = data[, idx, drop = FALSE],
            kinship   = K[[chrom]],
            genoprobs = snp_prob,
            addcovar  = covar,
            intcovar  = interactive_covariate,
            cores     = num_cores
        )

        map_tmp <- snpinfo_to_map(window_snps)
        tmp <- expand_snp_results(out_snps, map_tmp, window_snps)

        ret$lod_intcovar <- tmp$lod[, 1]
    }

    ret
}



# snpinfo to map
# direct copy from rqtl/qtl2
snpinfo_to_map <-
    function(snpinfo)
{
    uindex <- sort(unique(snpinfo$index))
    if(any(snpinfo$index < 1 | snpinfo$index > nrow(snpinfo)))
        stop("snpinfo$index values outside of range [1, ",
             nrow(snpinfo), "]")

    uchr <- unique(snpinfo$chr)
    chr <- factor(snpinfo$chr, levels=uchr)

    map <- split(snpinfo$pos, chr)
    snp <- split(snpinfo$snp, chr)
    index <- split(snpinfo$index, chr)
    for(i in seq_along(map)) {
        u <- unique(index[[i]])
        map[[i]] <- map[[i]][u]
        names(map[[i]]) <- snp[[i]][u]
    }

    names(map) <- uchr

    map
}


# expand snp association results according to snpinfo
# direct copy from rqtl/qtl2
expand_snp_results <- function(snp_results, map, snpinfo) {
    snpinfo <- split(snpinfo, factor(snpinfo$chr, unique(snpinfo$chr)))

    if(length(map) != length(snpinfo))
        stop("length(map) [", length(map), "] != length(snpinfo) [",
             length(snpinfo), "]")

    if(nrow(snp_results) != length(unlist(map)))
        stop("nrow(snp_results) [", nrow(snp_results),
             "] != length(unlist(map)) [",
             length(unlist(map)), "]")

    cnames <- rep(names(map), vapply(map, length, 0))
    lodindex <-
        split(seq_len(nrow(snp_results)), factor(cnames, unique(cnames)))

    result <- NULL
    for(i in seq(along=map)) {
        revindex <- rev_snp_index(snpinfo[[i]])

        map[[i]] <- snpinfo[[i]]$pos
        names(map[[i]]) <- snpinfo[[i]]$snp
        this_result <-
            unclass(snp_results)[lodindex[[i]],,drop=FALSE][revindex,,drop=FALSE]
        rownames(this_result) <- snpinfo[[i]]$snp

        result <- rbind(result, this_result)
    }

    list(lod=result,
         map=map)
}

# reverse index
# direct copy from rqtl/qtl2
rev_snp_index <- function(snpinfo) {
    index_spl <- split(seq_len(nrow(snpinfo)), snpinfo$index)
    revindex <- rep(seq_along(index_spl), vapply(index_spl, length, 1))
    revindex[unlist(index_spl)] <- revindex

    revindex
}
