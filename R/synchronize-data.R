#' Synchronize data and subset accordingly.
#'
#' @param dataset the dataset object
#'
#' @return list with 3 elements: `annots`, `samples` and `data`.
#'
#' @export
synchronize_data <- function(dataset) {
    # get the annotations from the dataset
    annotations <- get_annotations(dataset)
    annotation_ids <- annotations$annotation_id

    # get the samples from the dataset
    samples <- get_samples(dataset)
    sample_ids <- samples$sample_id

    # grab the data, row names are sample ids, column names are annotation ids
    data <- get_data(dataset)

    # see where all the sample ids are
    samples_only_in_data <- setdiff(rownames(data), sample_ids)
    samples_only_in_samples <- setdiff(sample_ids, rownames(data))
    sample_ids <- intersect(sample_ids, rownames(data))

    if (length(sample_ids) == 0) {
        stop("There are no samples in common")
    }

    # see where all the annotation ids are
    annotations_only_in_data <- setdiff(colnames(data), annotation_ids)
    annotations_only_in_annots <- setdiff(annotation_ids, colnames(data))
    annotation_ids <- intersect(annotation_ids, colnames(data))

    if (length(annotation_ids) == 0) {
        message("There are no annotations in common")
    }

    # sort the ids to make sure the come back in order
    annotation_ids <- sort(annotation_ids)
    sample_ids <- sort(sample_ids)

    # filter the annotations
    annotations <- annotations %>%
        dplyr::filter(.data$annotation_id %in% annotation_ids)

    # filter the data
    data <- data[sample_ids, annotation_ids, drop = FALSE]

    # filter the samples
    samples <- samples %>% dplyr::filter(.data$sample_id %in% sample_ids)

    list(
        annotations                  = annotations,
        #annotations_only_in_data     = annotations_only_in_data,
        #annotations_only_in_annots   = annotations_only_in_annots,
        data                         = data,
        samples                      = samples
        #samples_only_in_data         = samples_only_in_data,
        #samples_only_in_samples      = samples_only_in_samples
    )
}
