#' Validate the data in the dataset to make sure qtl2api can use it.
#'
#' @param dataset the dataset_id as a string identifier
#'
#' @export
validate_data <- function(dataset) {
    cat("STATUS  : Checking data\n")

    # get the dataset
    if (is.character(dataset)) {
        ds_orig <- get_dataset_by_id(dataset)
    } else {
        ds_orig <- dataset
    }

    ds <- synchronize_dataset(ds_orig)

    samples_name <- grep("^samples?$",
                         names(ds_orig),
                         value = TRUE)

    num_annot_samples <- nrow(ds_orig[[samples_name]])


    annotations_name <- grep("^annots?$|annotations?$",
                             names(ds_orig),
                             value = TRUE)

    annotations_orig <- ds_orig[[annotations_name]]

    # again allow some play in the column name
    annotations_id_field <- grep(
        "^annots?(\\.|_){0,1}ids?$|^annotations?(\\.|_){0,1}ids?$",
        colnames(annotations_orig),
        value = TRUE,
        ignore.case = TRUE
    )

    annotation_ids <- ds[[annotations_name]][[annotations_id_field]]

    if (is_phenotype(ds)) {
        annots_temp <-
            ds[[annotations_name]] %>%
            dplyr::filter(
                .data$omit == FALSE,
                .data$is_pheno == TRUE
            )
        annotation_ids <- annots_temp[[annotations_id_field]]
    }

    if ("data" %not in% names(ds)) {
        message("ERROR   : data not found")
        return()
    }

    if (is.matrix(ds$data)) {
        # check if the data is numeric
        if (!is.numeric(ds$data)) {
            message("ERROR   : data matrix is not numeric")
            return()
        }

        if (num_annot_samples != nrow(ds$data)) {
            cat("WARNING : number of samples (", num_annot_samples, ") != ",
                "number of data rows (", nrow(ds$data), ")\n")
        }

        if (NCOL(ds$data) != NROW(annotation_ids)) {
            cat('WARNING : number of annotations (', NROW(annotation_ids), ') != ',
                'number of data cols (', NCOL(ds$data), ")\n")

            x <- setdiff(annotation_ids, colnames(ds$data))
            y <- setdiff(colnames(ds$data), annotation_ids)

            if (length(x) > 0) {
                #cat("WARNING : annotations with no data:", paste(x, sep = ",", collapse = ","), "\n")
                cat("WARNING : # annotations with no data:", length(x), "\n")
            }

            if (length(y) > 0) {
                #cat("WARNING : data with no annotations:", paste(y, sep = ",", collapse = ","), "\n")
                cat("WARNING : # data with no annotations:", length(y), "\n")
            }
        }

        #perc_missing = (sum(is.na(ds_orig$data)) * 100) / prod(dim(ds_orig$data))
        #cat("STATUS  : Percentage of missing original data: ", perc_missing, "\n")

        perc_missing = (sum(is.na(ds$data)) * 100) / prod(dim(ds$data))
        if (perc_missing >= 10) {
            cat("WARNING : Percentage of missing synchronized data: ", perc_missing, "\n")
        } else {
            cat("STATUS  : Percentage of missing synchronized data: ", perc_missing, "\n")
        }

    } else if (is.list(ds$data)) {
        # TODO: this won't happen on a synchronized dataset
        # no list, just a matrix for data
        data.found <- FALSE
        if (any(c('rz','norm','log','transformed','raw') %in% tolower(names(ds$data)))) {
            data.found <- TRUE
        }

        if (!data.found) {
            message("ERROR   : 'rz','norm','log','transformed', OR 'raw' NOT found in data, must be ONE of them")
        }

        data_list <- get('data', ds)
        data_names <- ls(data_list)

        for (i in 1:length(data_names)) {
            cat('checking ', data_names[i])
            data_to_check <- paste0('data$', data_names[i])
            temp_data <- get(data_names[i], data_list)
            if (!is.numeric(temp_data)) {
                message(paste0(data_to_check, ' is not numeric'))
            }

            if (num_annot_samples != nrow(temp_data)) {
                cat('WARNING : number of samples (', num_annot_samples, ') != ',
                    'number of data["', data_names[i], '"] rows (', nrow(temp_data), ')\n')
            }

            if (NCOL(temp_data) != NROW(annot_ids)) {
                cat('WARNING : number of annotations (', NROW(annot_ids), ') != ',
                               'number of data["', data_names[i], '"] cols (', NCOL(temp_data), ')\n')

                x <- setdiff(annot_ids, colnames(temp_data))
                y <- setdiff(colnames(temp_data), annot_ids)

                if (length(x) > 0) {
                    #cat("WARNING : annotations with no data:", paste(x, sep = ",", collapse = ","), "\n")
                    cat("WARNING : # annotations with no data:", length(x), "\n")
                }

                if (length(y) > 0) {
                    #cat("WARNING : data with no annotations:", paste(y, sep = ",", collapse = ","), "\n")
                    cat("WARNING : # data with no annotations:", length(y), "\n")
                }
            }

            #perc_missing = (sum(is.na(ds_orig$data[data_names[i]])) * 100) / prod(dim(ds_orig$data[data_names[i]]))
            #cat("STATUS  : Percentage of missing original data: ", perc_missing, "\n")

            perc_missing = (sum(is.na(temp_data)) * 100) / prod(dim(temp_data))
            cat("STATUS  : Percentage of missing synchronized data: ", perc_missing, "\n")
        }
    }
}
