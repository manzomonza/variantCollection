#' Aggregate all detected DNA variants per sample
#'
#' @param watchdogDir 
#'
#' @return
#' @export
#'
#' @examples
variantCollection <- function(watchdogDir){
  if(dir.exists(watchdogDir)){
    #print("dir exists")
    files = list.files(path = watchdogDir, full.names = TRUE)
    # grep different filepaths
    
    cnv = grep("prep_cnv", files, value = TRUE)
    snv = grep("prep_snv", files, value = TRUE)
    filtered = grep("prep_filtered", files, value = TRUE)
    
    file_list = list(cnv = cnv,
                     snv = snv,
                     filtered = filtered)
    dirpath = rev(stringr::str_split(dirname(watchdogDir), pattern = "/", simplify = TRUE))[1]
    
    ## Workaround due to messy Info.csv setup
    if(grepl("Snvindel_watchdog", watchdogDir)){
      #print("Dir is from precision")
      infopath = paste0(dirname(watchdogDir), "/Info.csv")
      if(file.exists(infopath)){
        #print('infopath exists')
        info = precisionInfo(infopath)
      }else{
        #print('infopath does not exist')
        info = data.frame(workflowName = "Precision",
                          analysisName = dirpath,
                          exportDate = NA,
                          analysisDate = NA)
      }
    }else{
      file_list$info = grep("Info.csv", files, value = TRUE)
      info = readInfo(file_list$info)
    }
    
    info$dirpath = dirpath
    #print('dirpath')
    cnv = dplyr::bind_cols(readr::read_tsv(file_list$cnv), info) %>%
      dplyr::mutate(across(.cols = everything(), .fns = as.character))
    #print('cnv pass')
    snv = dplyr::bind_cols(readr::read_tsv(file_list$snv), info) %>%
      dplyr::mutate(across(.cols = everything(), .fns = as.character))
    filtered = dplyr::bind_cols(readr::read_tsv(file_list$filtered), info) %>%
      dplyr::mutate(across(.cols = everything(), .fns = as.character))
    return(list(snv = snv,
                cnv = cnv,
                filtered = filtered))
  }
}
