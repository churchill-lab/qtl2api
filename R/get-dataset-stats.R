#' Get all "dataset.*" statistics.
#'
#' This will return a named list of all datasets and some statistics.
#'
#' @return A named list of all the dataset objects.
#' @export
get_dataset_stats <- function() {
    datasets <- utils::apropos('^dataset\\.*', ignore.case = TRUE)
    ret <- c()

    for (d in datasets) {
        ds <- get(d)

        annots_field <- NA

        if (tolower(ds$datatype) == 'mrna') {
            annots_field <- grep("^annots?(\\.|_){1}mrnas?$",
                                 names(ds),
                                 value = TRUE)
        } else if(tolower(ds$datatype) == 'protein') {
            annots_field <- grep("^annots?(\\.|_){1}proteins?$",
                                 names(ds),
                                 value = TRUE)
        } else if(tolower(ds$datatype) == 'protein_uniprot') {
            annots_field <- grep("^annots?(\\.|_){1}proteins?(\\.|_){1}uniprots?$",
                                 names(ds),
                                 value = TRUE)
        } else if(tolower(ds$datatype) == 'phos') {
            annots_field <- grep("^annots?(\\.|_){1}phos?$",
                                 names(ds),
                                 value = TRUE)
        } else if(tolower(ds$datatype) == 'ptm') {
            annots_field <- grep("^annots?(\\.|_){1}ptm?$",
                                 names(ds),
                                 value = TRUE)
        } else if(tolower(ds$datatype) == 'peptide') {
            annots_field <- grep("^annots?(\\.|_){1}peptide?$",
                                 names(ds),
                                 value = TRUE)
        } else if(is_phenotype(ds) == 'phenotype') {
            annots_field <- grep("^annots?(\\.|_){1}pheno(type)?s?$",
                                 names(ds),
                                 value = TRUE)
        }

        annots_field_samples <- grep("^annots?(\\.|_){1}samples?$",
                                     names(ds),
                                     value = TRUE)

        display_name_field <- grep(
            "^display(\\.|_){1}name$",
            names(ds),
            ignore.case = TRUE,
            value = TRUE
        )

        display_name <- d

        if (length(display_name_field) != 0) {
            display_name <- ds[[display_name_field]]
        }

        temp <- list(id              = d,
                     display_name    = display_name,
                     datatype        = ds$datatype,
                     annotations     = annots_field,
                     num_annotations = NROW(ds[[annots_field]]),
                     samples         = annots_field_samples,
                     num_samples     = NROW(ds[[annots_field_samples]]))

        ret <- c(ret, list(temp))
    }

    ret
}
