#' Get the annotations from the dataset.
#'
#' @param dataset the dataset object
#'
#' @return the annotations tibble
#'
#' @export
get_annotations <- function(dataset) {
    # allow some play in the name annot, annots, annotation, annotations
    annotations_name <- grep("^annots?$|annotations?$",
                             names(dataset),
                             value = TRUE)

    if (length(annotations_name) == 0) {
        stop("Unable to find any annotations tibble in dataset")
    }

    annotations <-
        dataset[[annotations_name]] %>%
        janitor::clean_names()

    if (is_phenotype(dataset)) {
        annotations <-
            annotations %>%
            dplyr::filter(.data$omit == FALSE, .data$is_pheno == TRUE)
    }

    # again allow some play in the column name
    annotations_id_field <- grep(
        "^annots?(\\.|_){0,1}ids?$|^annotations?(\\.|_){0,1}ids?$",
        colnames(annotations),
        value = TRUE,
        ignore.case = TRUE
    )

    # make sure we pass back a tibble with annotation_id
    if (annotations_id_field != "annotation_id") {
        annotations$annotation_id <- annotations[[annotations_id_field]]
    }

    # make sure annotation_id is the first column (easier to view)
    annotations <-
        annotations %>%
        dplyr::relocate(annotation_id = .data$annotation_id)

    if (!is_phenotype(dataset)) {
        # make sure we are using base pairs
        annotations$start <- annotations$start %>% tidyr::replace_na(0)
        annotations$end <- annotations$end %>% tidyr::replace_na(0)

        annotations <- annotations %>%
            dplyr::mutate(
                start = as.integer(ifelse(
                    .data$start < 1000,
                    .data$start * 1000000,
                    .data$start
                ))
            ) %>%
            dplyr::mutate(
                end = as.integer(ifelse(
                    .data$end < 1000,
                    .data$end * 1000000,
                    .data$end
                ))
            )
    }

    return (annotations)
}
