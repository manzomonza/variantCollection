## VARIANT COLLECTION SOURCE

# snv_path = '/Users/manzo/Downloads/Sample_centric_SNV.tsv'
# cnv_path = '/Users/manzo/Downloads/Sample_centric_CNV.tsv'
# filtered_path = '/Users/manzo/Downloads/Sample_centric_FILTERED.tsv'

snv_path = '/mnt/NGS_Diagnostik/SampleCentricVariantCollection/Sample_centric_SNV.tsv'
cnv_path = '/mnt/NGS_Diagnostik/SampleCentricVariantCollection/Sample_centric_CNV.tsv'
filtered_path = '/mnt/NGS_Diagnostik/SampleCentricVariantCollection/Sample_centric_FILTERED.tsv'
variant_ints_path =  '/mnt/NGS_Diagnostik/SampleCentricVariantCollection/Variant_interpretations.tsv'


snv_parse <- function(snv_path){
  snv_table = readr::read_tsv(snv_path)
  snv_table = snv_table %>% 
    dplyr::filter(multiply_freq_by_100 != "multiply_freq_by_100")
  snv_table$percent_frequency = as.numeric(snv_table$percent_frequency)
  snv_table$percent_frequency = ifelse(snv_table$multiply_freq_by_100, snv_table$percent_frequency*100, snv_table$percent_frequency)
  snv_table = snv_table %>% dplyr::select(gene, locus, coding,
                                          one_AA, percent_frequency,
                                          analysisName, analysisDate,workflowName)
  snv_table = snv_table %>% dplyr::rename(one_AminoAcid_change = one_AA,
                                          gene_symbol = gene)
  snv_table = snv_table %>% dplyr::filter(analysisDate != 'exon')
  snv_table = dplyr::distinct(snv_table)
  return(snv_table)
}

cnv_parse <- function(cnv_path){
  cnv_table = readr::read_tsv(cnv_path)
  cnv_table = cnv_table %>%
    dplyr::rename( gene_symbol = gene) %>%
    dplyr::select(gene_symbol, locus, copy_number, contains("Percent"),analysisName, analysisDate, workflowName)
  
  cnv_table$fivePercent_conf = as.numeric(cnv_table$fivePercent_conf)
  cnv_table$ninetyfivePercent_conf = as.numeric(cnv_table$ninetyfivePercent_conf)
  
  cnv_table = cnv_table %>% dplyr::filter(!is.na(ninetyfivePercent_conf) & !is.na(fivePercent_conf))
  cnv_table = dplyr::distinct(cnv_table)
  return(cnv_table)
}

filtered_parse <- function(filtered_path){
  filtered_table = readr::read_tsv(filtered_path)
  filtered_table = filtered_table %>%
    dplyr::select(gene,locus, coding, contains("amino"), percent_frequency,analysisName,analysisDate, workflowName)
  filtered_table = filtered_table %>% dplyr::filter(!is.na(coding) & !is.na(amino_acid_change))
    filtered_table = dplyr::distinct(filtered_table)
  return(filtered_table)
}

variant_interpretations_parse <- function(variant_ints_path){
  variants = readr::read_tsv(variant_ints_path)
  variants = variants %>% dplyr::filter(!is.na(COSMIC_n_total) | 
                                        !is.na(COSMIC_n_tissue),
                                        !is.na(BIMI_variant) | 
                                        !is.na(Kommentar) | 
                                        !is.na(tsgInfo) |
                                        !is.na(cancerHotspot))
  variants$gene = ifelse(is.na(variants$gene), variants$genes, variants$gene)
  variants = variants %>% dplyr::select(gene, coding, one_AA, COSMIC_n_total, cancerHotspots, BIMI_variant, Kommentar )
  return(variants)
}