#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(markdown)
library(knitr)
library(shinymaterial)
library(shinyjs)


# Define UI for application that draws a histogram
ui <- material_page(
  useShinyjs(),
  
  tags$head(
    includeCSS("lib/custom_theme.css")
  ),
  # Application title
  title = "ProTN: integrative pipeline for the analysis of proteomics data from MS", 
  
  navbarPage(
    windowTitle = "ProTN",
    # Sidebar with a slider input for number of bins 
    # splitLayout(cellArgs = list(style='white-space: normal;'),
    # splitLayout(cellWidths = c("60%", "20%", "20%"), cellArgs = list(style='white-space: normal;'),
        # Show a plot of the generated distribution
        # mainPanel(
          # uiOutput("mark")
    tagList(
      div(style="display: flex; flex-direction: row; align-items: flex-end;",
      material_modal(
        modal_id = "execution",
        button_text = "Run the analysis",
              
        title = "ProTN",
        tagList(
          material_text_box('title_exp','Title of Analysis'),
          material_text_box('description_exp','Brief description'),
          material_radio_button('sw_analyzer', 'Software Analyzer', c('PD', 'MQ')),
          material_file_input('input_file', 'Select the INPUT file...'),
          material_file_input('pep_file', 'Select the PEP file...'),
          material_file_input('prot_file', 'Select the PROT file...'),
          material_text_box('signal_DEPs','Signal log2 expr thr', value = "inf"),
          material_text_box('FC_DEPs','Log2 FC thr',value = "0.75"),
          material_text_box('pvalue_DEPs','P.Value thr',value = "0.05"),
          material_radio_button('batch_corr', 'Batch effect correction', c('TRUE', 'FALSE'), selected = FALSE),
          material_text_box('prot_boxplot','Control Boxplot proteins'),
          material_file_input('design', 'Select a file with the design for the comparisons...'),
          material_radio_button('STRING', 'Execute PPI network STRINGdb', c('TRUE', 'FALSE'), selected = FALSE),
          material_radio_button('enrichR', 'Execute enrichment', c('TRUE', 'FALSE'), selected = FALSE),
          uiOutput("input_enrichment"),
          downloadButton("report", "Generate report"),
          plotOutput("n_plot")
        )
      ),
      material_button('contact', 'Contacts')
    ))
  ),
  
  includeMarkdown("README.md")
          
          # HTML(markdownToHTML(knit("README.md"), fragment.only=TRUE))
            # plotOutput("distPlot")
        # )
        # mainPanel(
        #   textInput('title_exp','Title of Analysis'),
        #   textAreaInput('description_exp','Brief description', rows = 6),
        #   radioButtons('sw_analyzer', 'Software Analyzer', c('PD', 'MQ'), inline = TRUE),
        #   fileInput('input_file', 'Select the INPUT file...'),
        #   fileInput('pep_file', 'Select the PEP file...'),
        #   fileInput('prot_file', 'Select the PROT file...'),
        #   textInput('signal_DEPs','Signal log2 expr thr', value = "inf"),
        #   textInput('FC_DEPs','Log2 FC thr',value = "0.75"),
        #   textInput('pvalue_DEPs','P.Value thr',value = "0.05")
        # ),
        # mainPanel(
        #   radioButtons('batch_corr', 'Batch effect correction', c('TRUE', 'FALSE'), inline = TRUE, selected = FALSE),
        #   textInput('prot_boxplot','Control Boxplot proteins'),
        #   fileInput('design', 'Select a file with the design for the comparisons...'),
        #   radioButtons('STRING', 'Execute PPI network STRINGdb', c('TRUE', 'FALSE'), inline = TRUE, selected = FALSE),
        #   radioButtons('enrichR', 'Execute enrichment', c('TRUE', 'FALSE'), inline = TRUE, selected = FALSE),
        #   uiOutput("input_enrichment"),
        #   downloadButton("report", "Generate report")
        # )
    # )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize=10*1024^3)
  
  #Visibility of the enrichment parameters based on the value of the Enrichment radiobutton
  output$input_enrichment<-renderUI({
    if (input$enrichR){
      tagList(
        material_text_box('pvalue_enrich','P.Value thr for enrichment',value = "0.05"),
        material_slider("os_enrich", "Overlap size thr for enrichment", 1, 30, step_size = 1, initial_value = 5),
        material_text_box('terms_enrich','Terms to search'),
        selectizeInput('DB_enrich', 'DB to analyse', choices = read_delim("lib/dbs_enrichR.txt", delim = "\n", col_names = FALSE)[[1]], 
                       selected = NULL, multiple = TRUE)
      )
    } 
  })
  
  #READ PARAMETERS
  output$report <- downloadHandler( 
    filename = "results.zip",
    content = function(file) {
      tryCatch({
        material_spinner_show(session, "n_plot")
        dirOutput_2=tempdir()
        currentTime = gsub(".*?([0-9]+).*?", "\\1", Sys.time())
        dirOutput_1=paste("/",currentTime,"/",sep = "")
        dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
        dirOutput_Server=paste(dirOutput_2,dirOutput_1,sep = "")
        message(dirOutput_Server)
        # tempReport <- file.path(tempdir(), "lib/pipeline_elaborate_PD_files.Rmd")
        # file.copy("pipeline_elaborate_PD_files.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(doc_title = input$title_exp,
                       description = input$description_exp,
                       readPD_files = if(input$sw_analyzer =="PD"){TRUE}else{FALSE},
                       readMQ_files = if(input$sw_analyzer =="MQ"){TRUE}else{FALSE},
                       file_input = input$input_file$datapath,
                       file_prot = input$prot_file$datapath,
                       file_pep = input$pep_file$datapath,
                       signal_thr = input$signal_DEPs,
                       fc_thr = input$FC_DEPs,
                       pval_thr = input$pvalue_DEPs,
                       batch_corr_exe = input$batch_corr,
                       contr_design = input$design$datapath,
                       prot_boxplot = input$prot_boxplot,
                       run_enrich = input$enrichR,
                       run_STRING = input$STRING,
                       pval_enrich_thr = input$pvalue_enrich,
                       overlap_size_enrich_thr = input$os_enrich,
                       enrich_filter_term = input$terms_enrich,
                       enrich_filter_DBs = input$DB_enrich,
                       dirOutput = dirOutput_Server)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render("lib/pipeline_elaborate_PD_files.Rmd", output_file = "report.html",
                          output_dir = dirOutput_Server,
                          params = params,
                          envir = new.env(parent = globalenv()))
                          
        oldwd<-getwd()
        setwd(dirOutput_Server)
        files2zip <- list.files("./", recursive = TRUE)
        #TODO
        zip(zipfile = file, files = files2zip, extra="-r")
        setwd(oldwd)
        material_spinner_hide(session, "n_plot")
        
      }, error = function(e) {
        showNotification(paste0("ERROR: ",e), type = "error")
      })
    }
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server, options = list(port = 8100, host = "127.0.0.1"))
