library(shiny)
library(DT)

source("SampleCentricTables_paths.R")

refreshtime = 3600000
### SHINY PART

ui <- fluidPage(
  titlePanel("Aggregation of detected variants, CNVs and filtered entries."),
  fluidRow(
    column(12,
           p("Data is derived from watchdog folder entries in the respective cases."))
  ),
  tabsetPanel(id = "tabs",
              tabPanel(value = "snv", title = "SNV",
                       DTOutput('snv')
              ),
              tabPanel(value = "cnv", title = "CNV",
                       DTOutput("cnv")
              ),
              tabPanel(value = "filtered", title = "Filtered",
                       DTOutput("filtered")
              ),
              tabPanel(value = "variants_int", title = "Variant Interpretations",
                       DTOutput("variants_int")
              )
  )
)
server <- function(input, output , session) {
  snvReader <- reactiveFileReader(refreshtime, session, snv_path, snv_parse)
  cnvReader <- reactiveFileReader(refreshtime, session, cnv_path, cnv_parse)
  filteredReader <- reactiveFileReader(refreshtime, session, filtered_path, filtered_parse)
  variant_int_Reader <- reactiveFileReader(refreshtime, session, variant_ints_path, variant_interpretations_parse)
  
  output$snv <- renderDT(snvReader(),
                         filter = "top",
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10,25,-1),
                                                          c(10,25,"All")))
  )
  output$cnv <- renderDT(cnvReader(),
                         filter = "top",
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10,25,-1),
                                                          c(10,25,"All")))
  )
  output$filtered <- renderDT(filteredReader(),
                              filter = "top",
                              extensions = 'Buttons',
                              options = list(dom = 'Blfrtip',
                                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                             lengthMenu = list(c(10,25,-1),
                                                               c(10,25,"All")))
  )
  output$variants_int <- renderDT(variant_int_Reader(),
                              filter = "top",
                              extensions = 'Buttons',
                              options = list(dom = 'Blfrtip',
                                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                             lengthMenu = list(c(10,25,-1),
                                                               c(10,25,"All")))
  )
}

shinyApp(ui, server)
