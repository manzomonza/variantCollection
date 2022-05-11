#  watchdog output based variant aggregation:
# Strategy becomes to extract as much as possible from watchdog folder, including metadata. Known exceptions are:
# Precision panels:
#   - will need to be checked one directory above watchdog
#   - maybe does not have info.csv --> Extract filename from path, rest metadata NA
# State: April 2022

library(magrittr)
library(optparse)


filteredFile <- "/home/ionadmin/ngs_variant_annotation/variantCollection/Sample_centric_FILTERED.tsv"
snvFile <- "/home/ionadmin/ngs_variant_annotation/variantCollection/Sample_centric_SNV.tsv"
cnvFile <- "/home/ionadmin/ngs_variant_annotation/variantCollection/Sample_centric_CNV.tsv"


### To be selected columns
filtered_cols <- c("amino_acid_change","analysisDate","analysisName","cds_region","cds_region_No",
                   "cnv_confidence", "coding", "copy_number", "dirpath",
                   "exon", "exportDate", "gene","location","locus","multiply_freq_by_100",
                   "percent_frequency","transcript","type","workflowName")
cnv_cols <- c("analysisDate","analysisName","chromosome","copy_number","dirpath","exportDate",
              "five","fivePercent_conf","gene","locus","ninetyfive","ninetyfivePercent_conf", 
              "workflowName")

snv_cols <- c("amino_acid_change","analysisDate","analysisName","cds_region","cds_region_No","clinvar_ready_AA",
              "coding","dirpath","exon","exportDate","gene","IR_clinvar","location","locus","multiply_freq_by_100",
              "one_AA","percent_frequency","three_AA","transcript","type","workflowName")



### OPT Parse
option_list = list(
  make_option(c("-d", "--watchdir"), type="character", default=NULL,
              help="watchdog directory", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

print(opt$watchdir)

## Function call
metavariants <- variantCollection(opt$watchdir)

## Select only columns shared by all panels (April 2022)
metavariants$snv <- metavariants$snv %>% dplyr::select(all_of(snv_cols))
metavariants$cnv <- metavariants$cnv %>% dplyr::select(all_of(cnv_cols))
metavariants$filtered <- metavariants$filtered %>% dplyr::select(all_of(filtered_cols))

# Write out -- SNV
writeVariants(tableIn = metavariants$snv,
              tableOut = snvFile)
# Write out --  CNV
writeVariants(tableIn = metavariants$cnv,
              tableOut = cnvFile)
# Write out --  FILTERED
writeVariants(tableIn = metavariants$filtered,
              tableOut = snvFiltered)
