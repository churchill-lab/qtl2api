' Summarize LOD scores for cis, distal, and trans associations.
#'
#' Given LOD scores at markers for a single gene and that gene's annotation,
#' this function extracts the maximum LOD score, and also determines the
#' LOD score for cis, distal, and trans associations.
#'
#' If no markers fall within the cis window, the nearest marker on the same
#' chromosome is used as a proxy for cis.
#'
#' @param lod_scores_mod A tibble with LOD results. Must contain columns: `id`, `chr`, `pos`, and `lod`.
#' @param annot_info A list or named vector with gene information. Must contain `chr` and `start`.
#' @param cis_window The cis region size in base pairs (default: 10 Mb).
#'
#' @return A named `list` with the following elements:
#' \itemize{
#'   \item max_lod - Maximum LOD score overall
#'   \item max_marker - Marker ID corresponding to max LOD
#'   \item cis_lod - Max LOD score within cis window
#'   \item cis_marker - Marker ID for cis LOD
#'   \item distal_lod - Max LOD score on same chromosome outside cis
#'   \item distal_marker - Marker ID for distal LOD
#'   \item trans_lod - Max LOD score on other chromosomes
#'   \item trans_marker - Marker ID for trans LOD
#' }
lod_summarize <- function(lod_scores_mod, annot_info, cis_window = 10e6) {

    # Find the index of the maximum LOD score across all markers
    max_idx <- which.max(lod_scores_mod$lod)

    # Extract the maximum LOD score and the corresponding marker ID
    max_lod <- lod_scores_mod$lod[max_idx]
    max_marker <- lod_scores_mod$id[max_idx]

    # Initialize cis, distal, and trans values as NA
    cis_lod <- NA
    cis_marker <- NA
    distal_lod <- NA
    distal_marker <- NA
    trans_lod <- NA
    trans_marker <- NA

    # Continue only if annotation has necessary information
    if (!is.null(annot_info) && ('chr' %in% names(annot_info)) && ('start' %in% names(annot_info))) {

        # Define the safe boundaries of the cis window, making sure the start doesn't go below 0
        cis_lower <- max(0, annot_info$start - cis_window)
        cis_upper <- annot_info$start + cis_window

        # Logical vector: markers on the same chromosome
        is_same_chr <- lod_scores_mod$chr == annot_info$chr

        # Logical vector: markers within cis window on the same chromosome
        is_cis <- is_same_chr & (lod_scores_mod$pos >= cis_lower & lod_scores_mod$pos <= cis_upper)

        # Logical vector: markers on the same chromosome but outside the cis window
        is_distal <- is_same_chr & !is_cis

        # Logical vector: markers on other chromosomes
        is_trans <- !is_same_chr

        ## ---- CIS LOGIC ----
        if (any(is_cis)) {
            # If there are cis markers, get the one with the highest LOD score
            cis_idx <- which.max(lod_scores_mod$lod[is_cis])
            cis_lod <- lod_scores_mod$lod[is_cis][cis_idx]
            cis_marker <- lod_scores_mod$id[is_cis][cis_idx]
        } else if (any(is_same_chr)) {
            # Fallback: if no markers in cis window, find nearest marker on same chromosome
            dists <- abs(lod_scores_mod$pos[is_same_chr] - annot_info$start)
            nearest_idx <- which.min(dists)
            same_chr_rows <- which(is_same_chr)
            nearest_row <- same_chr_rows[nearest_idx]
            cis_lod <- lod_scores_mod$lod[nearest_row]
            cis_marker <- lod_scores_mod$id[nearest_row]
        }

        ## ---- DISTAL LOGIC ----
        if (any(is_distal)) {
            distal_idx <- which.max(lod_scores_mod$lod[is_distal])
            distal_lod <- lod_scores_mod$lod[is_distal][distal_idx]
            distal_marker <- lod_scores_mod$id[is_distal][distal_idx]
        }

        ## ---- TRANS LOGIC ----
        if (any(is_trans)) {
            trans_idx <- which.max(lod_scores_mod$lod[is_trans])
            trans_lod <- lod_scores_mod$lod[is_trans][trans_idx]
            trans_marker <- lod_scores_mod$id[is_trans][trans_idx]
        }
    }

    # Return all results in a named list
    return(list(
        max_lod = max_lod,
        max_marker = max_marker,
        cis_lod = cis_lod,
        cis_marker = cis_marker,
        distal_lod = distal_lod,
        distal_marker = distal_marker,
        trans_lod = trans_lod,
        trans_marker = trans_marker
    ))
}


#' Perform a LOD scan.
#'
#' @param dataset The dataset object.
#' @param id The unique id in the dataset.
#' @param intcovar The interactive covariate.
#' @param cores Number of cores to use (0 = ALL).
#' @param filter_threshold If set, qtl2::find_peaks is used.
#' @param filter_peak_drop If set, qtl2::find_peaks is used.
#' @param filter_thresholdX If set, qtl2::find_peaks is used.
#' @param filter_peak_dropX If set, qtl2::find_peaks is used.
#' @param scan1_output If `TRUE`, original `qtl2::scan1` data is included.
#' @param cis_window cis defined window
#'
#' @return a `list` with the following elements:
#' \itemize{
#'   \item lod_peaks - `tibble` of LOD peaks
#'   \item lod_scores - `tibble` with the following columns: id, chr, pos, lod
#'   \item scan1 - `qtl2::scan1` output
#' }
#'
#' @importFrom rlang .data
#' @export
calc_lod_scores <- function(dataset, id, intcovar = NULL, cores = 0,
                            filter_threshold = 6.0, filter_peak_drop = Inf,
                            filter_thresholdX = NULL, filter_peak_dropX = NULL,
                            scan1_output = FALSE, cis_window = 10e6) {
    # make sure samples and annotations are available
    ds <- synchronize_dataset(dataset)

    # check if id exists
    if (!any(id == colnames(ds$data))) {
        stop(sprintf("Cannot find id '%s' in dataset", id))
    }

    annot_info <- NULL

    if (tolower(ds$datatype) == "mrna") {
        annot_info <-
            ds$annot_mrna %>%
            dplyr::filter(.data$gene_id == id)
    } else if (tolower(ds$datatype) == "protein") {
        annot_info <-
            ds$annot_protein %>%
            dplyr::filter(.data$protein_id == id)
    } else if (tolower(ds$datatype) == "protein_uniprot") {
        annot_info <-
            ds$annot_protein_uniprot %>%
            dplyr::filter(.data$uniprot_id == id)
    } else if (tolower(ds$datatype) == "phos") {
        annot_info <-
            ds$annot_phos %>%
            dplyr::filter(.data$phos_id == id)
    } else if (tolower(ds$datatype) == "ptm") {
        annot_info <-
            ds$annot_ptm %>%
            dplyr::filter(.data$ptm_id == id)
    } else if (tolower(ds$datatype) == "peptide") {
        annot_info <-
            ds$annot_peptide %>%
            dplyr::filter(.data$peptide_id == id)
    } else if (tolower(ds$datatype) == "phenotype") {
        annot_info <-
            ds$annot_phenotype %>%
            dplyr::filter(.data$data_name == id)
    }

    if (is.null(annot_info) || nrow(annot_info) == 0) {
        stop(sprintf("Cannot find annotation for id '%s' in dataset", id))
    } else {
        # for now, just take the first one
        annot_info <-
            annot_info %>%
            dplyr::slice(1) %>%
            as.list()
    }

    # make sure num_cores is appropriate
    num_cores <- nvl_int(cores, 0)

    # get the covar information
    covar_information <- get_covar_matrix(ds, id)
    covar_matrix <- covar_information$covar_matrix
    covar_formula <- covar_information$covar_formula

    # set the interactive.covariate, to be used in scan1
    interactive_covariate <- NULL

    if (!is.null(intcovar)) {
        if (!any(intcovar == ds$covar_info$sample_column)) {
            stop(sprintf("intcovar '%s' not found in covar_info", intcovar))
        }

        if (is.null(covar_matrix)) {
            stop(sprintf("no covar_matrix, but intcovar '%s' specified", intcovar))
        }

        # grabbing all the columns from covar (covar.matrix) that
        # match, i.e., "batch" will match "batch2", "BATCH3", etc
        interactive_covariate <-
            covar_matrix[, which(grepl(intcovar, colnames(covar_matrix), ignore.case = T))]
    }

    # perform the scan using QTL2,
    # - addcovar should always be ALL covars
    # - intcovar should be just the interactive covariate column
    lod_scores <- qtl2::scan1(
        genoprobs = genoprobs,
        kinship   = K,
        pheno     = ds$data[, id, drop = FALSE],
        addcovar  = covar_matrix,
        intcovar  = interactive_covariate,
        cores     = num_cores,
        reml      = TRUE
    )

    # utilize qtl2::find_peaks
    lod_peaks <- qtl2::find_peaks(
        lod_scores,
        map,
        threshold  = filter_threshold,
        peakdrop   = filter_peak_drop,
        thresholdX = filter_thresholdX,
        peakdropX  = filter_peak_dropX
    )

    markers_cleaned <-
        markers %>%
        dplyr::filter(!is.na(.data$pos)) %>%
        janitor::clean_names()

    # construct a 2 dimensional array of data with id, chr, pos, lod as columns
    # we perform a left join here to make sure that the number of elements match
    # also convert from type scan1 to numeric
    lod_scores_mod <-
        dplyr::inner_join(
            tibble::as_tibble(lod_scores, rownames = "marker_id"),
            markers_cleaned,
            by = "marker_id"
        ) %>%
        dplyr::select(
            id  = .data$marker_id,
            chr = .data$chr,
            pos = .data$pos,
            lod = id
        ) %>%
        dplyr::mutate_at(
            c("lod"),
            as.numeric
        )

    lod_peaks <-
        dplyr::inner_join(
            lod_peaks,
            markers_cleaned,
            by = c("chr", "pos")
        ) %>%
        dplyr::select(
            id  = .data$marker_id,
            chr = .data$chr,
            pos = .data$pos,
            lod = .data$lod
        ) %>%
        dplyr::mutate_at(
            c("lod"),
            as.numeric
        ) %>%
        tibble::as_tibble()

    # convert to bp
    if (all(lod_scores_mod$pos < 1000)) {
        lod_scores_mod$pos <- as.integer(lod_scores_mod$pos * 1000000)
    }

    if (all(lod_peaks$pos < 1000)) {
        lod_peaks$pos <- as.integer(lod_peaks$pos * 1000000)
    }

    lod_ext <- lod_summarize(lod_scores_mod, annot_info, cis_window = cis_window)


    ret <- list(
        lod_peaks  = lod_peaks,
        lod_scores = lod_scores_mod,
        annot_info = annot_info,
        lod_ext    = lod_ext
    )

    ret$scan1 <- if (scan1_output) lod_scores else NULL

    attr(ret, 'covar_formula') <- covar_formula

    ret
}
