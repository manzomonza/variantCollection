## VARIANT COLLECTION SOURCE

snv_table = readr::read_tsv('/mnt/NGS_Diagnostik/SampleCentricVariantCollection/Sample_centric_SNV.tsv')
cnv_table = readr::read_tsv('/mnt/NGS_Diagnostik/SampleCentricVariantCollection/Sample_centric_CNV.tsv')  %>% dplyr::distinct()
filtered_table = readr::read_tsv('/mnt/NGS_Diagnostik/SampleCentricVariantCollection/Sample_centric_FILTERED.tsv')  %>%
  dplyr::distinct() 