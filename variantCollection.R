#/usr/bin/R
### Scope
### Use NGSannotation methods to read in files
### and export metadata at once

library(magrittr)
library(NGSannotation)
library(dplyr)

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
      
      print(cnv_filepath)
      print(snv_filepath)
      cnv = readIn(cnv_filepath)
      snv = read_rename_select_precision_snv(snv_filepath)
      output_tables <- make_output_tables_precision(snv_indel = snv, cnv = cnv)
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
  return()
}


aggregator(opt$file)







