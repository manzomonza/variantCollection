#/usr/bin/R
### Scope
### Use NGSannotation methods to read in files
### and export metadata at once

library(magrittr)
library(NGSannotation)
library(tidyverse)
### and export metadata at once
library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

print(opt$file)
### 
aggregator <- function(filepath){
  if(nrow(readIn(filepath) > 0)){
    if(grepl("Snvindel.tsv", filepath, ignore.case = FALSE)){
      cnv_filepath =  paste0(dirname(filepath),"/Cnv.tsv")
      snv_filepath =  paste0(dirname(filepath),"/Snvindel.tsv")
      info_csv = paste0(dirname(filepath),"/Info.csv")
      if(file.exists(info_csv)){
        metadata_foi = precisionInfo(info_csv)
      }else{
        filename = rev(stringr::str_split(dirname(filepath), pattern = "/", simplify = TRUE))[1]
        metadata_foi = tibble(analysisName =  metadata_foi,
                      analysisDate = NA,
                      exportDate = NA,
                      workflowName = NA)   
      }
      print(cnv_filepath)
      print(snv_filepath)
      cnv = readIn(cnv_filepath)
      snv = read_rename_select_precision_snv(snv_filepath)
      output_tables <- make_output_tables_precision(snv_indel = snv, cnv = cnv)
      output_tables = lapply(output_tables, function(x) dplyr::bind_cols(metadata_foi, x))
      
    }else{
      ir_output = read_rename_select(filepath)
      output_tables = make_output_tables(ir_output)
      metadata_foi = metadataCollection(filepath)
      output_tables = lapply(output_tables, function(x) dplyr::bind_cols(metadata_foi, x))
      }
  return(output_tables)
    }
}
metadataCollection <- function(filepath){
  metadat = readr::read_tsv(filepath, col_names = F) %>%
    dplyr::filter(grepl("##", X1)) %>%
    dplyr::mutate(X1 = gsub("##", "", X1))
  colnames(metadat) = "IR_Workflow_Metainformation"
  metadat = metadat %>% tidyr::separate(col = IR_Workflow_Metainformation, into = c("parameter", "value"), sep = "=") %>%
    tidyr::pivot_wider(names_from = parameter, values_from = value)
  return(metadat)
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

metavariants = aggregator(opt$file)


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


