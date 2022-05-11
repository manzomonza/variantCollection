#' Read metadata from Precision
#'
#' @param infopath 
#'
#' @return Parsed metadata as data.frame
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
