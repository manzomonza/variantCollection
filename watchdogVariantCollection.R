# # watchdog output based variant aggregation:
# # Strategy becomes to extract as much as possible from watchdog folder, including metadata. Known exceptions are:
# # Precision panels:
# #   - will need to be checked one directory above watchdog
# #   - maybe does not have info.csv --> Extract filename from path, rest metadata NA
# # State: April 2022
library(magrittr)
library(optparse)


### OPT Parse
option_list = list(
  make_option(c("-d", "--watchdir"), type="character", default=NULL,
              help="watchdog directory", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

print(opt$watchdir)


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


precisionInfo <- function(filepath){
  info = readr::read_csv(filepath, col_names = TRUE) %>%
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


variantCollection <- function(watchdogDir){
  if(dir.exists(watchdogDir)){
    #print("dir exists")
    files <- list.files(path = watchdogDir, full.names = TRUE)
    # grep different filepaths
    file_list <- list(cnv = grep("prep_cnv", files, value = TRUE),
                      snv = grep("prep_snv", files, value = TRUE),
                      filtered = grep("prep_filtered", files, value = TRUE))
    dirpath <- rev(stringr::str_split(dirname(watchdogDir), pattern = "/", simplify = TRUE))[1]
    
    
    if(grepl("Snvindel_watchdog", watchdogDir)){
      #print("Dir is from precision")
      infopath = paste0(dirname(watchdogDir), "/Info.csv")
      if(file.exists(infopath)){
        #print('infopath exists')
        info = precisionInfo(infopath)
      }else{
        #print('infopath does not exist')
        info = data.frame(workflowName = "Precision")
      }
    }else{
      file_list$info = grep("Info.csv", files, value = TRUE)
      info = readInfo(file_list$info)
    }
    
    info$dirpath = dirpath
    #print('dirpath')
    cnv = dplyr::bind_cols(readr::read_tsv(file_list$cnv), info) %>% dplyr::mutate(across(.cols = everything(), .fns = as.character))
    #print('cnv pass')
    snv = dplyr::bind_cols(readr::read_tsv(file_list$snv), info) %>% dplyr::mutate(across(.cols = everything(), .fns = as.character))
    filtered = dplyr::bind_cols(readr::read_tsv(file_list$filtered), info) %>% dplyr::mutate(across(.cols = everything(), .fns = as.character))
    return(list(snv = snv,
                cnv = cnv,
                filtered = filtered))
  }
}

snv_cols <- c("amino_acid_change","analysisDate","analysisName","cds_region","cds_region_No","clinvar_ready_AA",
              "coding","dirpath","exon","exportDate","gene","IR_clinvar","location","locus","multiply_freq_by_100",
              "one_AA","percent_frequency","three_AA","transcript","type","workflowName")

cnv_cols <- c("analysisDate","analysisName","chromosome","copy_number","dirpath","exportDate",
              "five","fivePercent_conf","gene","locus","ninetyfive","ninetyfivePercent_conf", 
              "workflowName")

filtered_cols <- c("amino_acid_change","analysisDate","analysisName","cds_region","cds_region_No",
                   "cnv_confidence", "coding", "copy_number", "dirpath",
                   "exon", "exportDate", "gene","location","locus","multiply_freq_by_100",
                   "percent_frequency","transcript","type","workflowName")

## Function call
metavariants <- variantCollection(opt$watchdir)

## Select only columns shared by all panels (April 2022)
metavariants$snv <- metavariants$snv %>% dplyr::select(all_of(snv_cols))
metavariants$cnv <- metavariants$cnv %>% dplyr::select(all_of(cnv_cols))
metavariants$filtered <- metavariants$filtered %>% dplyr::select(all_of(filtered_cols))


if(nrow(metavariants$snv) > 0){
  if(file.exists("Sample_centric_SNV.tsv")){
    readr::write_tsv(metavariants$snv, "Sample_centric_SNV.tsv", append = TRUE, col_names = TRUE)
  }else{
    readr::write_tsv(metavariants$snv, "Sample_centric_SNV.tsv", append = FALSE, col_names = TRUE)
  }
}

if(nrow(metavariants$cnv) > 0){
  if(file.exists("Sample_centric_CNV.tsv")){
    readr::write_tsv(metavariants$cnv, "Sample_centric_CNV.tsv", append = TRUE, col_names = TRUE)
  }else{
    readr::write_tsv(metavariants$cnv, "Sample_centric_CNV.tsv", append = FALSE, col_names = TRUE)
  }
}

if(nrow(metavariants$filtered) > 0){
  if(file.exists("Sample_centric_FILTERED.tsv")){
    readr::write_tsv(metavariants$filtered, "Sample_centric_FILTERED.tsv", append = TRUE, col_names = TRUE)
  }else{
    readr::write_tsv(metavariants$filtered, "Sample_centric_FILTERED.tsv", append = FALSE, col_names = TRUE)
  }
}




