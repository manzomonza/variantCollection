library(shiny)
library(DT)

source("SampleCentricTables_paths.R")

## SNV
snv_table = snv_table %>% 
  dplyr::filter(multiply_freq_by_100 != "multiply_freq_by_100")
snv_table$percent_frequency = as.numeric(snv_table$percent_frequency)
snv_table$percent_frequency = ifelse(snv_table$multiply_freq_by_100, snv_table$percent_frequency*100, snv_table$percent_frequency)
snv_table = snv_table %>% dplyr::select(gene, locus, coding,
                                        one_AA, percent_frequency,
                                        analysisName, workflowName, analysisDate)
snv_table = snv_table %>% dplyr::rename(one_AminoAcid_change = one_AA,
                                        gene_symbol = gene)
snv_table = snv_table %>% dplyr::distinct()

## CNV
cnv_table = cnv_table %>%
  dplyr::rename( gene_symbol = gene) %>%
  dplyr::select(gene_symbol, locus, copy_number, contains("Percent"), workflowName, analysisDate)

cnv_table$fivePercent_conf = as.numeric(cnv_table$fivePercent_conf)
cnv_table$ninetyfivePercent_conf = as.numeric(cnv_table$ninetyfivePercent_conf)

cnv_table = cnv_table %>% dplyr::filter(!is.na(ninetyfivePercent_conf) & !is.na(fivePercent_conf))


## FILTERED
filtered_table = filtered_table %>%
  dplyr::select(gene,locus, coding, contains("amino"), percent_frequency,
                analysisName, workflowName) %>%
  dplyr::filter(!is.na(coding) & !is.na(amino_acid_change))

### SHINY PART

ui <- fluidPage(
  tabsetPanel(id = "tabs",
              tabPanel(value = "snv", title = "SNV",
                       DTOutput('snv')
              ),
              tabPanel(value = "cnv", title = "CNV",
                       DTOutput("cnv")
              ),
              tabPanel(value = "filtered", title = "Filtered",
                       DTOutput("filtered")
              )
  )
)
server <- function(input, output) {
  output$snv <- renderDT(snv_table,
                         filter = "top",
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10,25,-1),
                                                          c(10,25,"All")))
  )
  output$cnv <- renderDT(cnv_table,
                         filter = "top",
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10,25,-1),
                                                          c(10,25,"All")))
  )
  output$filtered <- renderDT(filtered_table,
                              filter = "top",
                              extensions = 'Buttons',
                              options = list(dom = 'Blfrtip',
                                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                             lengthMenu = list(c(10,25,-1),
                                                               c(10,25,"All")))
  )
  
}

shinyApp(ui, server)
