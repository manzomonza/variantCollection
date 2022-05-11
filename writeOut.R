#' Write out results
#'
#' @param tableIn is a dataframe containing all
#' @param tableOut is character string designating table destination
#'
#' @return
#' @export
#'
#' @examples
writeVariants <- function(tableIn, tableOut){
  if(nrow(tableIn) > 0){
    if(file.exists(tableOut)){
      readr::write_tsv(tableIn, tableOut, append = TRUE, col_names = TRUE)
    }else{
      readr::write_tsv(tableIn, tableOut, append = FALSE, col_names = TRUE)
    }
  }
}
