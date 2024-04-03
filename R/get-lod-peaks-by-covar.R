#' Get the LOD peaks for a dataset.
#'
#' Used internally to the package.
#'
#' By default, this will get all LOD peaks stored in the `dataset$lod_peaks`
#' `list`.  If `intcovar` is specified, only the peaks for that are returned.
#'
#' @param ds The dataset object to get the LOD peaks for.
#' @param intcovar NULL for 'additive' or a value representing the interactive
#' covariate.
#'
#' @return A a `tibble` with the following columns depending upon
#' `dataset$datatype`:
#' \itemize{
#'   \item mRNA - `marker,chr,pos,gene_id,symbol,gene_chrom,middle,lod`
#'   \item protein - `marker,chr,pos,protein_id,gene_id,symbol,gene_chrom,middle,lod`
#'   \item phenotype - `marker,chr,pos,data_name,short_name,description,lod`
#' }
#'
#' If the allele effects are stored, values `A-H` will also be included.
#' @export
get_lod_peaks_by_covar <- function(ds, intcovar = NULL) {
    # these functions should not be synchronized
    if (valid(ds$is_synchronized)) {
        stop("dataset should not be synchronized")
    }

    annots_field_peaks <- grep("^lod(\\.|_){1}peaks?$",
                               names(ds),
                               value = TRUE)

    if (length(annots_field_peaks) == 0) {
        return(NULL)
    } else if (is.null(ds[[annots_field_peaks]])) {
        return(NULL)
    }

    annots_field_covar <- grep("^covar(\\.|_){1}info$",
                               names(ds),
                               value = TRUE)

    peaks <- NULL

    if ((is.null(intcovar)) || (intcovar == 'additive')) {
        peaks <- ds[[annots_field_peaks]]$additive
    } else {
        # find the covar and get the name of the lod peaks
        covar_info <- ds[[annots_field_covar]] %>% janitor::clean_names()

        if (any(intcovar == covar_info$sample_column)) {
            n <- covar_info[covar_info$sample_column == intcovar, ]
            peaks <- ds[[annots_field_peaks]][[n$lod_peaks]]
        }
    }

    if (invalid(peaks)) {
        stop(sprintf("No peaks found for intcovar '%s' in lod.peaks", intcovar))
    }

    # this is a little extra work because we are trying to be nice for users
    # who separate with '.' or '_'
    peaks %<>% janitor::clean_names()

    markers_cleaned <-
        markers %>%
        dplyr::filter(!is.na(.data$pos)) %>%
        janitor::clean_names()

    # convert from Mbp to bp
    if (all(markers_cleaned$pos < 1000)) {
        markers_cleaned$pos <- as.integer(markers_cleaned$pos * 1000000)
    }

    if (tolower(ds$datatype) == "mrna") {
        annots_field <- grep("^annots?(\\.|_){1}mrnas?$",
                             names(ds),
                             value = TRUE)

        annots <- ds[[annots_field]] %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names()

        if (all(annots$start < 1000)) {
            annots$start <- as.integer(annots$start * 1000000)
        }

        if (all(annots$end < 1000)) {
            annots$end <- as.integer(annots$end * 1000000)
        }

        ret <- annots %>%
            dplyr::inner_join(
                peaks,
                by = "gene_id"
            ) %>%
            dplyr::select(
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                gene_chr   = .data$chr,
                start      = .data$start,
                end        = .data$end,
                marker_id  = .data$marker_id,
                lod        = .data$lod
            ) %>%
            dplyr::mutate(
                gene_pos = as.integer(round((.data$start + .data$end) / 2))
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                lod        = .data$lod,
                marker_id  = .data$marker_id,
                chr        = .data$chr,
                pos        = .data$pos,
                gene_chr   = .data$gene_chr,
                gene_pos   = .data$gene_pos,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
            ) %>%
            dplyr::arrange(
                .data$chr,
                .data$pos
            )

        # now add A-H for additive if they exist
        if (all(tolower(LETTERS[1:8]) %in% colnames(peaks))) {
            ret <- ret %>%
                dplyr::inner_join(
                    peaks,
                    by = c("gene_id", "marker_id", "lod")
                )
        }
    } else if (tolower(ds$datatype) == "protein") {
        annots_field <- grep("^annots?(\\.|_){1}proteins?$",
                             names(ds),
                             value = TRUE)

        annots <- ds[[annots_field]] %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names()

        if (all(annots$start < 1000)) {
            annots$start <- as.integer(annots$start * 1000000)
        }

        if (all(annots$end < 1000)) {
            annots$end <- as.integer(annots$end * 1000000)
        }

        ret <- annots %>%
            dplyr::inner_join(
                peaks,
                by = "protein_id"
            ) %>%
            dplyr::select(
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                gene_chr   = .data$chr,
                start      = .data$start,
                end        = .data$end,
                marker_id  = .data$marker_id,
                lod        = .data$lod
            ) %>%
            dplyr::mutate(
                gene_pos = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                lod        = .data$lod,
                marker_id  = .data$marker_id,
                chr        = .data$chr,
                pos        = .data$pos,
                gene_chr   = .data$gene_chr,
                gene_pos   = .data$gene_pos,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                protein_id = .data$protein_id
            ) %>%
            dplyr::arrange(
                .data$chr,
                .data$pos
            )

        # now add A-H for additive if they exist
        if (all(tolower(LETTERS[1:8]) %in% colnames(peaks))) {
            ret <- ret %>%
                dplyr::inner_join(
                    peaks,
                    by = c("protein_id", "marker_id", "lod")
                )
        }
    } else if (tolower(ds$datatype) == "protein_uniprot") {
        annots_field <- grep("^annots?(\\.|_){1}proteins?(\\.|_){1}uniprots?$",
                             names(ds),
                             value = TRUE)

        annots <- ds[[annots_field]] %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names()

        #annots$start <- annots$start %>% tidyr::replace_na(0)
        #annots$end <- annots$end %>% tidyr::replace_na(0)

        if (all(annots$start < 1000)) {
            annots$start <- as.integer(annots$start * 1000000)
        }

        if (all(annots$end < 1000)) {
            annots$end <- as.integer(annots$end * 1000000)
        }

        ret <- annots %>%
            dplyr::inner_join(
                peaks,
                by = "uniprot_id"
            ) %>%
            dplyr::select(
                uniprot_id = .data$uniprot_id,
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                gene_chr   = .data$chr,
                start      = .data$start,
                end        = .data$end,
                marker_id  = .data$marker_id,
                lod        = .data$lod
            ) %>%
            dplyr::mutate(
                gene_pos = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                lod        = .data$lod,
                marker_id  = .data$marker_id,
                chr        = .data$chr,
                pos        = .data$pos,
                gene_chr   = .data$gene_chr,
                gene_pos   = .data$gene_pos,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                protein_id = .data$protein_id,
                uniprot_id = .data$uniprot_id
            ) %>%
            dplyr::arrange(
                .data$chr,
                .data$pos
            )

        # now add A-H for additive if they exist
        if (all(tolower(LETTERS[1:8]) %in% colnames(peaks))) {
            ret <- ret %>%
                dplyr::inner_join(
                    peaks,
                    by = c("uniprot_id", "marker_id", "lod")
                )
        }
    } else if (tolower(ds$datatype) == "phos") {
        annots_field <- grep("^annots?(\\.|_){1}phos$",
                             names(ds),
                             value = TRUE)

        annots <- ds[[annots_field]] %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names()

        #annots$start <- annots$start %>% tidyr::replace_na(0)
        #annots$end <- annots$end %>% tidyr::replace_na(0)

        if (all(annots$start < 1000)) {
            annots$start <- as.integer(annots$start * 1000000)
        }

        if (all(annots$end < 1000)) {
            annots$end <- as.integer(annots$end * 1000000)
        }

        ret <- annots %>%
            dplyr::inner_join(
                peaks,
                by = "phos_id"
            ) %>%
            dplyr::select(
                phos_id    = .data$phos_id,
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                gene_chr   = .data$chr,
                start      = .data$start,
                end        = .data$end,
                marker_id  = .data$marker_id,
                lod        = .data$lod
            ) %>%
            dplyr::mutate(
                gene_pos = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                lod        = .data$lod,
                marker_id  = .data$marker_id,
                chr        = .data$chr,
                pos        = .data$pos,
                gene_chr   = .data$gene_chr,
                gene_pos   = .data$gene_pos,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                protein_id = .data$protein_id,
                phos_id.   = .data$phos_id
            ) %>%
            dplyr::arrange(
                .data$chr,
                .data$pos
            )

        # now add A-H for additive if they exist
        if (all(tolower(LETTERS[1:8]) %in% colnames(peaks))) {
            ret <- ret %>%
                dplyr::inner_join(
                    peaks,
                    by = c("phos_id", "marker_id", "lod")
                )
        }
    } else if (tolower(ds$datatype) == "peptide") {
        annots_field <- grep("^annots?(\\.|_){1}peptide$",
                             names(ds),
                             value = TRUE)

        annots <- ds[[annots_field]] %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names()

        #annots$start <- annots$start %>% tidyr::replace_na(0)
        #annots$end <- annots$end %>% tidyr::replace_na(0)

        if (all(annots$start < 1000)) {
            annots$start <- as.integer(annots$start * 1000000)
        }

        if (all(annots$end < 1000)) {
            annots$end <- as.integer(annots$end * 1000000)
        }

        ret <- annots %>%
            dplyr::inner_join(
                peaks,
                by = "peptide_id"
            ) %>%
            dplyr::select(
                peptide_id = .data$peptide_id,
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                uniprot_id = .data$uniprot_id,
                gene_chr   = .data$chr,
                start      = .data$start,
                end        = .data$end,
                marker_id  = .data$marker_id,
                lod        = .data$lod
            ) %>%
            dplyr::mutate(
                gene_pos = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                lod        = .data$lod,
                marker_id  = .data$marker_id,
                chr        = .data$chr,
                pos        = .data$pos,
                gene_chr   = .data$gene_chr,
                gene_pos   = .data$gene_pos,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                protein_id = .data$protein_id,
                peptide_id = .data$peptide_id,
                uniprot_id = .data$uniprot_id
            ) %>%
            dplyr::arrange(
                .data$chr,
                .data$pos
            )

        # now add A-H for additive if they exist
        if (all(tolower(LETTERS[1:8]) %in% colnames(peaks))) {
            ret <- ret %>%
                dplyr::inner_join(
                    peaks,
                    by = c("peptide_id", "marker_id", "lod")
                )
        }
    } else if (tolower(ds$datatype) == "ptm") {
        annots_field <- grep("^annots?(\\.|_){1}ptm$",
                             names(ds),
                             value = TRUE)

        annots <- ds[[annots_field]] %>%
            dplyr::filter(!is.na(.data$chr) & !is.na(.data$start) & !is.na(.data$end)) %>%
            janitor::clean_names()

        #annots$start <- annots$start %>% tidyr::replace_na(0)
        #annots$end <- annots$end %>% tidyr::replace_na(0)

        if (all(annots$start < 1000)) {
            annots$start <- as.integer(annots$start * 1000000)
        }

        if (all(annots$end < 1000)) {
            annots$end <- as.integer(annots$end * 1000000)
        }

        ret <- annots %>%
            dplyr::inner_join(
                peaks,
                by = "ptm_id"
            ) %>%
            dplyr::select(
                ptm_id     = .data$ptm_id,
                peptide_id = .data$peptide_id,
                protein_id = .data$protein_id,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                uniprot_id = .data$uniprot_id,
                gene_chr   = .data$chr,
                start      = .data$start,
                end        = .data$end,
                marker_id  = .data$marker_id,
                lod        = .data$lod
            ) %>%
            dplyr::mutate(
                gene_pos = round((.data$start + .data$end) / 2)
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                lod        = .data$lod,
                marker_id  = .data$marker_id,
                chr        = .data$chr,
                pos        = .data$pos,
                gene_chr   = .data$gene_chr,
                gene_pos   = .data$gene_pos,
                gene_id    = .data$gene_id,
                symbol     = .data$symbol,
                protein_id = .data$protein_id,
                peptide_id = .data$peptide_id,
                ptm_id     = .data$ptm_id,
                uniprot_id = .data$uniprot_id
            ) %>%
            dplyr::arrange(
                .data$chr,
                .data$pos
            )

        # now add A-H for additive if they exist
        if (all(tolower(LETTERS[1:8]) %in% colnames(peaks))) {
            ret <- ret %>%
                dplyr::inner_join(
                    peaks,
                    by = c("ptm_id", "marker_id", "lod")
                )
        }
    } else if (is_phenotype(ds)) {
        annots_field <- grep("^annots?(\\.|_){1}pheno(type)?s?$",
                             names(ds),
                             value = TRUE)

        annots <- ds[[annots_field]] %>% janitor::clean_names()

        ret <- annots %>%
            dplyr::inner_join(
                peaks,
                by = "data_name"
            ) %>%
            dplyr::select(
                data_name   = .data$data_name,
                short_name  = .data$short_name,
                description = .data$description,
                marker_id   = .data$marker_id,
                lod         = .data$lod
            ) %>%
            dplyr::inner_join(
                markers_cleaned,
                by = "marker_id"
            ) %>%
            dplyr::select(
                marker_id   = .data$marker_id,
                chr         = .data$chr,
                pos         = .data$pos,
                data_name   = .data$data_name,
                short_name  = .data$short_name,
                description = .data$description,
                lod         = .data$lod
            ) %>%
            dplyr::arrange(
                .data$chr,
                .data$pos
            )

        # now add A-H for additive if they exist
        if (all(tolower(LETTERS[1:8]) %in% colnames(peaks))) {
            ret <- ret %>%
                dplyr::inner_join(
                    peaks,
                    by = c("data_name", "marker_id", "lod")
                )
        }
    } else {
        stop(sprintf("dataset has an invalid datatype '%s'", ds$datatype))
    }

    ret
}


