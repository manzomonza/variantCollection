#/usr/bin/R
### Scope
### Use NGSannotation methods to read in files
### and export metadata at once

library(magrittr)
library(NGSannotation)
library(tidyverse)
### and export metadata at once
library(optparse)
library(GenomicRanges)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

print(opt$file)
########## 


#### OVERWRITE NGS ANNOTATION function do not need snv exon annotation
read_rename_select_precision_snv <- function(filepath){
  frequency_col_mult100 = readIn(filepath) %>%
    mult100()
  snvindelx <- readIn(filepath) %>%
    addMissingCols() %>%
    renameCol() %>%
    snvParse() %>%
    dplyr::ungroup() %>%
    dplyr::select(-transcript) %>%
    dplyr::filter(grepl("PRESENT", call)) %>%
    dplyr::mutate(multiply_freq_by_100 = frequency_col_mult100) %>%
    dplyr::mutate(percent_frequency = as.numeric(percent_frequency))
  
  snvindelx <- exonAnnot(snvindelx) %>%
    dplyr::select(type, gene, coding, amino_acid_change, percent_frequency, location,
                  locus, transcript, exon, copy_number, cnv_confidence, IR_clinvar, multiply_freq_by_100) 
  return(snvindelx)
}

######## 

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
        metadata_foi = tibble::tibble(analysisName =  filename,
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



exon.gr <- ""
metavariants = aggregator("/Volumes/GoogleDrive/.shortcut-targets-by-id/1yuFiN1dlcUgo1_ELdNVXegTfB61oDv8G/Patientendaten/2022/W0501-W0550/W0510_PrecisionDNA/Snvindel.tsv")

precisionInfo()



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


