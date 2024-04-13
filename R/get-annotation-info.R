#' Get the annotation_info from the dataset.
#'
#' @param dataset the dataset object
#'
#' @return the annotation_info tibble
#'
#' @export
get_annotation_info <- function(dataset) {
    function(dataset) {
        # allow some play in the name annot.info, annots_info, annotation_info, etc
        annotation_info_name <-
            grep("^annots?(\\.|_){0,1}infos?$|^annotations?(\\.|_){0,1}infos?$",
                 names(dataset),
                 value = TRUE
            )

        if (length(annotation_info_name) == 0) {
            stop("Unable to find any annotation_info tibble in dataset")
        }

        annotation_info <-
            dataset[[annotation_info_name]] %>%
            janitor::clean_names()

        return (annotation_info)
    }
}
