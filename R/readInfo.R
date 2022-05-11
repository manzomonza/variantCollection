#' Read metadata from Info.csv
#'
#' @param infopath 
#'
#' @return
#' @export
#'
#' @examples
readInfo <- function(infopath){
  info = readr::read_csv(infopath) %>%
    tidyr::separate(col = 1,
                    into = c("No.", "metadata"), sep = "^\\d{1,} ") %>%
    tidyr::separate(col = "metadata",
                    into = c("metadata", "value"), sep = "=") %>%
    dplyr::select(-No.) %>%
    dplyr::mutate(value = stringr::str_remove_all(string = value, pattern = "NA| ")) %>%
    tidyr::pivot_wider(names_from = metadata,
                       values_from = value)
  if('chromosome' %in% colnames(info)){
    info = info %>% dplyr::rename(Chr = chromosome)
  }
  return(info)
}

#' Read metadata from Precision Info.csv
#'
#' @param infopath character string of path to Precision Info.csv
#'
#' @return
#' @export
#'
#' @examples
precisionInfo <- function(infopath){
  info = readr::read_csv(infopath, col_names = TRUE) %>%
    tidyr::separate(col = `Software Version Details`,
                    into = c("parameter", "value"),
                    sep = ",")
  info = info %>%
    as.data.frame() %>%
    dplyr::filter(grepl("Sample Name|^Name$|Start Date|Completion Date",parameter, ignore.case = TRUE)) %>%
    dplyr::filter(value != "<NA>" & value != "Assay Name" & value != "MT216761" & !grepl("gnxs--0036", value)) %>%
    tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
    dplyr::rename(analysisName = `Sample Name`,
                  workflowName = Name,
                  analysisDate = `Start Date`,
                  exportDate = `Completion Date`)
  return(info)
}

