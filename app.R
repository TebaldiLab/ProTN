################################################################################
# ProTN: an integrative pipeline for complete analysis of proteomics           # 
# data from mass spectrometry                                                  #
# Laboratory of RNA and Disease Data Science, University of Trento             #
# Developer: Gabriele Tomè                                                     #
# Issue at: https://github.com/TebaldiLab/ProTN/issues                         #
# PI: Dr. Toma Tebaldi, PhD                                                    #
################################################################################
list.of.packages <- c("shiny","tidyverse","markdown","knitr","shinydashboard",
                      "shinydashboardPlus","shinymaterial","shinyjs","magrittr",
                      "dplyr","stringr","shinyBS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) suppressMessages(suppressWarnings({install.packages(new.packages, dependencies = T)}))

#Load library
library(shiny)
library(tidyverse)
library(markdown)
library(knitr)
library(shinydashboard)
library(shinydashboardPlus)
library(shinymaterial)
library(shinyjs)
library(magrittr)
library(dplyr)
library(stringr)
library(shinyBS)

#Javascript code to disable click when is running the app
jsCode <- "shinyjs.pageDisable = function(params){
              $('body').css('pointer-events', params);
            };"

# Define UI for application that draws a histogram
ui <- tagList(
  dashboardPage(
    skin = "black",
    title = "ProTN",
    # Application title
    header=dashboardHeader(
      title=tags$a(href="#shiny-tab-info",
                   tags$img(id="logo_protn"))
    ),
    #Sidebar menu
    sidebar=dashboardSidebar(
      sidebarMenu(
        id="tabs",
        menuItem("Info ProTN", tabName = "info", icon = icon("info-circle", lib="font-awesome")),
        menuItem("Run ProTN", tabName = "analysis_protn", icon = icon("rocket", "fa-regular")),
        menuItem("Info PhosProTN", tabName = "info_phos", icon = icon("info-circle", lib="font-awesome")),
        menuItem("Run PhosProTN", tabName = "analysis_phosprotn", icon = icon("rocket", "fa-solid")),
        menuItem("Contacts", tabName = "contacts", icon = icon("comment", lib="font-awesome"))
      )
    ),
    #Main body web page
    body=dashboardBody(
      #Load associated file, CSS, JS, logo
      useShinyjs(),
      tags$meta(charset = "UTF-8"),
      tags$head(tags$link(rel="shortcut icon", href="images/logo_black.ico")),
      tags$head(tags$script(src="https://kit.fontawesome.com/5d5f342cf8.js")),
      tags$link(rel="stylesheet", href="https://fonts.googleapis.com/css?family=El+Messiri"),
      includeCSS("www/css/custom_theme.css"),
      includeCSS("www/css/materialize.css"),
      includeScript("www/js/materialize.js"),
      extendShinyjs(text = jsCode, functions = c("pageDisable")),
      
      #Busy panel when app is running
      conditionalPanel(
        condition = "$(\'html\').hasClass(\'shiny-busy\')",
        tags$div(class = "loader"),
        tags$div(class = "prevent_click")
      ),    
      tabItems(
        #INFO tab ProTN
        tabItem(
          tabName = "info",
          includeHTML("www/README.html"),
          bsModal("modal_Design_Example", 
                  "Design file Example", 
                  "preview_design",
                  size = "large",
                  DT::dataTableOutput("design_df")),
          bsModal("modal_Input_Example", 
                  "Input file Example", 
                  "preview_input",
                  size = "large",
                  DT::dataTableOutput("input_df")),
          bsModal("modal_peptide_Example", 
                  "Peptide file Example", 
                  "preview_peptide",
                  size = "large",
                  DT::dataTableOutput("pep_df")),
          bsModal("modal_protein_Example", 
                  "Protein file Example", 
                  "preview_protein",
                  size = "large",
                  DT::dataTableOutput("prot_df"))
        ),
        #Execution tab of ProTN
        tabItem(
          tabName = "analysis_protn",
          tagList(
            fluidRow(
              column(
                width = 4,
                fluidRow(
                  textInput("title_exp", "Title of the analysis"),
                  uiOutput("help1")
                ),
                fluidRow(
                  textAreaInput("description_exp", "Brief description", rows = 6),
                  uiOutput("help2")
                ),
                fluidRow(
                  radioButtons("sw_analyzer", "Software Analyzer", 
                               c("ProteomeDiscoverer", "MaxQuant", "TMT_PD"), inline = TRUE),
                  uiOutput("help3")
                )
              ),
              column(
                width = 4,
                fluidRow(
                  fileInput("input_file", "Select the SAMPLE_ANNOTATION file..."),
                  uiOutput("help4")
                ),
                fluidRow(
                  fileInput("pep_file", "Select the PEPTIDES file..."),
                  uiOutput("help5")
                ),
                fluidRow(
                  fileInput("prot_file", "Select the PROTEINS file..."),
                  uiOutput("help6")
                )
              ),
              column(
                width = 4,
                fluidRow(
                  fileInput("design", "Select a file with the design for the comparisons..."),
                  uiOutput("help12")
                ),
                radioButtons("custom_param", "Use custom parameter", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
                fillRow(
                  actionButton("preview_output_ProTN", "Generate report"),
                  actionButton("case_study_ProTN", "Case Study Example")
                )
              )
            ),
            uiOutput("input_params"),
            bsModal("modal_preview_case_study_ProTN", 
                    "Example case study", 
                    "case_study_ProTN",
                    size = "large",
                    downloadButton('download_report_CS', 'Download complete results (ZIP)'), uiOutput("preview_results_CS")),
            bsModal("modal_preview_output_ProTN", 
                    "Preview report of ProTN", 
                    "preview_output_ProTN",
                    size = "large",
                    downloadButton('download_report', 'Download complete results (ZIP)'), uiOutput("preview_results"))
          )
        ),
        #Info tab PhosProTN
        tabItem(
          tabName = "info_phos",
          includeHTML("www/README_phos.html"),
          bsModal("modal_Design_Example_phos", 
                  "Design file Example", 
                  "preview_design_phos",
                  size = "large",
                  DT::dataTableOutput("design_df_phos")),
          bsModal("modal_Input_Example_phos", 
                  "Input file Example", 
                  "preview_input_phos",
                  size = "large",
                  DT::dataTableOutput("input_df_phos")),
          bsModal("modal_peptide_Example_phos", 
                  "Peptide file Example", 
                  "preview_peptide_phos",
                  size = "large",
                  DT::dataTableOutput("pep_df_phos")),
          bsModal("modal_protein_Example_phos", 
                  "Protein file Example", 
                  "preview_protein_phos",
                  size = "large",
                  DT::dataTableOutput("prot_df_phos"))
        ),
        #Execution tab PhosProTN
        tabItem(
          tabName = "analysis_phosprotn",
          tagList(
            fluidRow(
              column(
                width = 4,
                fluidRow(
                  textInput("title_exp_phos", "Title of the analysis"),
                  uiOutput("help1_P")
                ),
                fluidRow(
                  textAreaInput("description_exp_phos", "Brief description", rows = 6),
                  uiOutput("help2_P")
                ),
                fluidRow(
                  radioButtons("sw_analyzer_phos", "Software Analyzer", c("ProteomeDiscoverer", "MaxQuant"), inline = TRUE),
                  uiOutput("help3_P")
                ),
                fluidRow(
                  fileInput("design_phos", "Select a file with the design for the comparisons..."),
                  uiOutput("help12_P")
                )
              ),
              column(
                width = 4,
                uiOutput("input_filePROT_phos")
              ),
              column(
                width = 4,
                uiOutput("input_filePHOS_phos"),
                actionButton("report_phos", "Generate report")
              )
            ),
            uiOutput("input_params_phos"),
            bsModal("modal_preview_output_PhosProTN", 
                    "Preview report of PhosProTN", 
                    "report_phos",
                    size = "large",
                    downloadButton('download_report_PhosProTN', 'Download complete results (ZIP)'), uiOutput("preview_results_phos"))
          )
        ),
        #Contact tab
        tabItem(
          tabName = "contacts",
          tagList(
            includeHTML("www/contacts.html")
          )
        )
      )
    )
  ),
  #Footer with lab, name and university
  tags$footer(
    tags$div(
      class="footer-row",
      tags$div(class="footer-col", style="text-align: left;",tags$p("This website is free and open to all users and there is no login requirement.")),
      tags$div(class="footer-col", style="text-align: right;",
               tags$div(
                 tags$img(id="logo_RDDS_footer",
                          src="data:image/svg+xml;base64,PHN2ZyBpZD0iTGF5ZXJfMSIgZGF0YS1uYW1lPSJMYXllciAxIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCA4MDUuMzMgODAxLjY3Ij48ZGVmcz48c3R5bGU+LmNscy0xe2ZpbGw6IzA4OGY4Zjt9LmNscy0ye2ZpbGw6bm9uZTtzdHJva2U6IzA4OGY4ZjtzdHJva2UtbWl0ZXJsaW1pdDoxMDtzdHJva2Utd2lkdGg6NXB4O308L3N0eWxlPjwvZGVmcz48ZyBpZD0iWWluWWFuZyI+PHBhdGggY2xhc3M9ImNscy0xIiBkPSJNNjM3LjE5LDQ2MS40N0E1LjIyLDUuMjIsMCwwLDAsNjQyLDQ2M2EyNy42LDI3LjYsMCwwLDAsMjMuMzctMTkuNDMsNS4xMiw1LjEyLDAsMCwwLS43My00Ljc5bC00Ny01Ni42MmEyMDUuNTUsMjA1LjU1LDAsMCwxLTI4LjY3LDIxLjE2WiIgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTI2LjgzIC0yMS42NykiLz48cGF0aCBjbGFzcz0iY2xzLTEiIGQ9Ik01NDMsNDFhMzc3LjU4LDM3Ny41OCwwLDAsMSwzOC40MSwxMS43OSw0MDMuMjQsNDAzLjI0LDAsMCwwLTU3LjcyLTE5LjFDMzA5LjM1LTE5Ljc5LDkyLjI1LDExMC42NSwzOC44MSwzMjVzNzcsNDMxLjQ0LDI5MS4zNSw0ODQuODlhNDAzLjY3LDQwMy42NywwLDAsMCw2Ni42NCwxMC44LDM3NS40NCwzNzUuNDQsMCwwLDEtNDYuMDktOC40OEMyNjMuNjIsNzkwLjQ5LDIwMS44OCw3MjIsMjAwLjgyLDYzNC4yN0gxMzMuNjRjLTYuNzMsMC0xMCwuNDUtMTMuMzctMi4zNnMtMi45My04LjQtMi45My0xMy4zNGMwLTguODEtLjQ2LTExLjExLDMuMzUtMTUsMy0zLDEwLjcxLTIuMDksMTctMiw2LjEzLjA5LDE0LjY4LDAsMTkuOTIsMCw2LjE0LDAsMTMtLjEzLDE5LjE3LDBzMTAuNCwwLDE2LjYzLDBoNS4xYTUuMzksNS4zOSwwLDAsMSwuNzEtLjA1LDYuNSw2LjUsMCwwLDEsLjc5LjA2aC4xOGE1Ljg0LDUuODQsMCwwLDEsLjczLS4wNSw2Ljg0LDYuODQsMCwwLDEsLjgxLjA2aDBhMTk0Ljg3LDE5NC44NywwLDAsMSw0LjUzLTI0LjkzYzExLjg3LTQ3LjYyLDQ4LTkyLjg1LDk1LjIxLTEyMi4wOUwyOTksNDUwLjE5bC01LjY5LTkuODVMMjg1LDQyNmwtNy4wOC0xMi4yN0wyNzAuNTUsNDAxYy0yLjYtNC41MS01LjUxLTguOTUtNy43OS0xMy42NGE1LjcyLDUuNzIsMCwwLDEtLjctMy4zMWMuMDUtMS4xLjE5LTUsMy43MS03LjQsMS4xMy0uNzYsMy43OC0yLjI5LDQuNC0yLjY1bDYtMy40NWMyLjU1LTEuNDcsNC44OS0zLjczLDkuMzMtMi44NSw2LDEuMiw4LjI0LDguMDUsMTEsMTIuOGwxMS4xMywxOS4yOCw0LjE0LDcuMTcsNC4yMiw3LjMyYzEuOTIsMy4zMSw3LjIxLDEyLjMzLDkuMDUsMTUuNjcuOSwxLjYxLDEuMzIsMi4xNCwxLjg0LDMuMTkuMTQuMjcuMjguNDIuNzUsMS4yOS42LDEuMTIuNzQsMS4xMiwxLjE0LDJsLjA5LjE4Yy40NS43NiwxLDEuODIsMS40MywyLjQ2LDM2Ljk1LTE2LjUsNzcuMzQtMjIuMjQsMTE2LjYxLTEyLjQ1QzQ5Ni41Niw0MzksNTQ3LDQyOSw1ODguOTIsNDAzLjM2bC01Ni40Ny02OGE0LjU2LDQuNTYsMCwwLDEtMS0zLjIyQTEyMiwxMjIsMCwwLDEsNDQxLDM0NS43OGMtNjUuNTEtMTUuNjItMTA2LjM3LTgyLjEyLTkwLjYzLTE0Ny42MWExMjIuMzgsMTIyLjM4LDAsMCwxLDE0OS43Mi04OS44OGM2NC4wNiwxNi42MSwxMDMuNTQsODIsODguNDEsMTQ2LjM5YTEyMS45MywxMjEuOTMsMCwwLDEtMzAuNjgsNTYuNjQsNC4yNiw0LjI2LDAsMCwxLDIsMS4zM2w1Ny43OCw2OS41NWEyMTAuMjgsMjEwLjI4LDAsMCwwLDY1LjA5LTEwNi44N0M3MDMsMTk0LDY0Niw2Ni42NSw1NDMsNDFaIiB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMjYuODMgLTIxLjY3KSIvPjxwYXRoIGNsYXNzPSJjbHMtMSIgZD0iTTU4MS40Miw1Mi43NUM3NjAuNzksMTI3LjU4LDg2My41LDMyNC4xOSw4MTUuMDUsNTE4LjU0Yy00Ny45LDE5Mi4xLTIyNy4yNCwzMTYuOC00MTguMjUsMzAyLjE1LDE5MC4zOCwyMi44OCwzNzMuOTQtMTA0LDQyMy4wOC0zMDEuMUM4NjkuNTksMzIwLjIxLDc2NC4xMiwxMTkuNjgsNTgxLjQyLDUyLjc1WiIgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTI2LjgzIC0yMS42NykiLz48L2c+PHBhdGggY2xhc3M9ImNscy0xIiBkPSJNNDY2LjM3LDgxMC4xOUEzNi40MiwzNi40MiwwLDAsMSw0NTIsODA3LjMyYTQzLjEyLDQzLjEyLDAsMCwxLTEzLjQ1LTkuMzgsMTAzLjE2LDEwMy4xNiwwLDAsMS0xMy43Mi0xNy42N2wtMjEuNjYtMzQuNTItOC41NiwwLC40NSw2Mi41LTI4LC4wNi0uOC0xMDkuNjhxLS4wOC05LjMxLTIuNzktMTVUMzUxLjA5LDY3OEgzNDlsMC02LjIxLDU5LjA4LS4xMnE1Ljc3LDAsMTMuMzguNjlhODcuNDUsODcuNDUsMCwwLDEsMTUuMzMsMi44Nyw1NS43Nyw1NS43NywwLDAsMSwxNC4zOSw2LjI4LDMyLjYxLDMyLjYxLDAsMCwxLDEwLjc4LDEwLjg1cTQuMTEsNi43MSw0LjE5LDE2LjY1LjA3LDExLTQuOSwxNy45MWEzNC42LDM0LjYsMCwwLDEtMTIuODcsMTAuODlBODUuNjYsODUuNjYsMCwwLDEsNDMxLjY0LDc0NGwyMi43NCwzNi4xN3E4LjI0LDEzLjIzLDE1LjExLDE3LjU2dDExLjYsNS4xNWwwLDUuMTdhMjMuMjYsMjMuMjYsMCwwLDEtNi4wOSwxLjQ2QTYxLjA4LDYxLjA4LDAsMCwxLDQ2Ni4zNyw4MTAuMTlabS03MS44MS03Mi43LDEyLjg0LDBBMzcuNzUsMzcuNzUsMCwwLDAsNDE3LjY2LDczNmEyNi4yMSwyNi4yMSwwLDAsMCw5LjI4LTQuNzgsMjMuNTgsMjMuNTgsMCwwLDAsNi42OC04LjkxLDMyLjYsMzIuNiwwLDAsMCwyLjQ3LTEzLjY3LDMzLjM0LDMzLjM0LDAsMCwwLTIuNTYtMTMuNTUsMjQuMjEsMjQuMjEsMCwwLDAtNi43LTguODgsMjUuODcsMjUuODcsMCwwLDAtOS4yNS00Ljg1LDM2LjgzLDM2LjgzLDAsMCwwLTEwLjE3LTEuNDJjLTEuMTUsMC0yLjgyLjA3LTUsLjIxYTYxLjMyLDYxLjMyLDAsMCwwLTguMjMsMS4yNloiIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0yNi44MyAtMjEuNjcpIi8+PHBhdGggY2xhc3M9ImNscy0xIiBkPSJNNTI5Ljc4LDc5NS4yOGwtNDYtMTAwLjIxcS0zLjkxLTguNTEtOC42MS0xMi42OXQtMTMuMDUtLjU0bC0xLjg1LjgxLTIuNi01LjY3LDUwLjQyLTIyLjA1cTM4LjE4LTE2LjcxLDY0LjkyLTguNjV0MzkuOTMsMzYuOHE4LjU5LDE4LjcsNi44NSwzNi41dC0xNCwzMi44cS0xMi4yNSwxNS0zNS42MSwyNS4yMlptMzYuMzgtMjVhNDMuNTcsNDMuNTcsMCwwLDAsMjAtMTYuNjdxNy0xMSw3LjA4LTI1LjhUNTg1LjQ4LDY5NnEtNy42My0xNi42My0xOC45My0yNi40NUE1MC44OCw1MC44OCwwLDAsMCw1NDIsNjU3LjM0YTQ1LjA3LDQ1LjA3LDAsMCwwLTI2LjU4LDMuNDVxLTQuMDgsMS43OS03LDMuMjdhODcuNjUsODcuNjUsMCwwLDAtNy44Miw0Ljc4bDQ5LjEzLDEwN3E0LjQzLTEsNi44Ny0xLjc2YTM3LjEzLDM3LjEzLDAsMCwwLDQuNDktMS42MloiIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0yNi44MyAtMjEuNjcpIi8+PHBhdGggY2xhc3M9ImNscy0xIiBkPSJNNjkxLjI3LDcwNi41OGwtODUuNjQtNzFxLTcuMjYtNi0xMy4zNS03LjgxdC0xMiw1LjA1TDU3OSw2MzQuMzVsLTQuODUtNEw2MDkuNzUsNTg5cTI3LTMxLjMxLDU0LjU1LTM1LjM5dDUyLjE2LDE2LjI4cTE2LDEzLjI3LDIyLjI5LDMwLjA4dDEuOTEsMzUuNTdxLTQuMzgsMTguNzYtMjAuODcsMzcuOTJabTIxLjY3LTM4QTQyLjU0LDQyLjU0LDAsMCwwLDcyMy41LDY0NXExLjQzLTEyLjktNS0yNi4zMnQtMjEtMjUuNDdxLTE0LjIyLTExLjc5LTI4LjcxLTE1Ljg2dC0yNy40NC0uNjFBNDQuMTMsNDQuMTMsMCwwLDAsNjE5LDU5MS4xNnEtMi44OCwzLjM1LTQuODEsNS45MWE4NS4zNyw4NS4zNywwLDAsMC00LjkyLDcuNjRsOTEuNDUsNzUuODJjMi4zNS0xLjg3LDQuMTUtMy4zOCw1LjQxLTQuNTFhMzguNTUsMzguNTUsMCwwLDAsMy4zMS0zLjM3WiIgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTI2LjgzIC0yMS42NykiLz48cGF0aCBjbGFzcz0iY2xzLTEiIGQ9Ik04MDguNzYsNTI2YTg2LjM0LDg2LjM0LDAsMCwxLTYuNzEsMTUuNzMsNjYuODQsNjYuODQsMCwwLDEtOS4zNCwxMy4xNiwzMy4wNywzMy4wNywwLDAsMS0xMSw4LjE3LDE2LDE2LDAsMCwxLTExLjQxLjc4bC0yMy42LTcuMjUsMS43Ni01LjZxMjIsMy4yMiwzMy45Mi0zLjc2dDE3LjIyLTIzLjgxYTM4LjE1LDM4LjE1LDAsMCwwLDItMTQsMTkuNjUsMTkuNjUsMCwwLDAtMy43OC0xMC44MywxOC43OSwxOC43OSwwLDAsMC05LjgyLTYuNDIsMjEsMjEsMCwwLDAtMTYuMTUsMS4wOSw3OS4xMSw3OS4xMSwwLDAsMC0xNi4zNiwxMS41NmwtMjIuNjcsMjBxLTEyLjM0LDExLjE1LTIyLjU5LDE0LjE1YTM2LjMyLDM2LjMyLDAsMCwxLTIxLjMzLS4zOXEtMTUuOTMtNC44OS0yMS40LTE5Ljg2dDEtMzUuNTlhODYuNjQsODYuNjQsMCwwLDEsNi43MS0xNS43Myw2Ni44NCw2Ni44NCwwLDAsMSw5LjM0LTEzLjE2LDM0LDM0LDAsMCwxLDEwLjg3LTguMiwxNiwxNiwwLDAsMSwxMS41Mi0uNzVsMjMuNiw3LjI0LTEuNzYsNS42MXEtMjItMy0zMy44OSwzLjY2VDY3OC4yOCw0ODMuOGEzMS4yMiwzMS4yMiwwLDAsMC0xLjUyLDEzLjI2QTE5LjIxLDE5LjIxLDAsMCwwLDY4MSw1MDcuMjVhMTkuNDIsMTkuNDIsMCwwLDAsOS41LDYsMTguNzMsMTguNzMsMCwwLDAsMTUuNTItMS41QTk5LjA3LDk5LjA3LDAsMCwwLDcyMS42Myw1MDBsMjIuNDgtMjAuMTFxMTIuMzUtMTEuMTQsMjIuMS0xNS4wN3QyMC44NC0uNTNxMTcuMzYsNS4zMywyMy4zNCwyMS40NFQ4MDguNzYsNTI2WiIgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTI2LjgzIC0yMS42NykiLz48cGF0aCBjbGFzcz0iY2xzLTEiIGQ9Ik00MjUuMjMsNDE1LjFoMzIuNjhhMCwwLDAsMCwxLDAsMHY4NC40NmE5LDksMCwwLDEtOSw5SDQzNC4yM2E5LDksMCwwLDEtOS05VjQxNS4xYTAsMCwwLDAsMSwwLDBaIiB0cmFuc2Zvcm09InRyYW5zbGF0ZSgyLjk3IC00OC40NCkgcm90YXRlKDMuNTkpIi8+PHBhdGggY2xhc3M9ImNscy0xIiBkPSJNMjQ5LjcyLDQ5My41M0gyODIuNGEwLDAsMCwwLDEsMCwwVjU3OGE5LDksMCwwLDEtOSw5SDI1OC43MmE5LDksMCwwLDEtOS05VjQ5My41M2EwLDAsMCwwLDEsMCwwWiIgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTM1My45MSA0MTQuMzcpIHJvdGF0ZSgtNTMuODIpIi8+PHBhdGggY2xhc3M9ImNscy0yIiBkPSJNMzQ2LjkyLDI3NC43OGM0LjUtNy4zNSw2LjgxLTksOC4wNS04LjU3LDEsLjM0LDEuMzksMi4xNiwyLjg2LDIuMzNhMy40OCwzLjQ4LDAsMCwwLDIuNi0xLjNjNC00LjA4LDkuODUtNS44MiwxNC41NC05LjA5LDEyLjgyLTguOTIsMzAuNjUtMS4zOSwzNC41Ni05LjYxLDEuMTctMi40NiwwLTQuMDgsMi4wNy03LjUzYTE1LjIxLDE1LjIxLDAsMCwxLDguNTgtNi41YzMuMDUtLjgxLDguMjktMS4wOCwxMS4xNywxLjgyLDEsMSwuOTMsMS42NywyLjYsMy4zOCwxLjM4LDEuNDEsNC4yNiwzLjg1LDUuNDUsMy4xMSwxLS42LS4zMi0yLjY0LDEtNC45M2E0LjU0LDQuNTQsMCwwLDEsMy4zOC0yLjM0YzEuNDQtLjA1LDEuODQsMS4wOCwyLjg2LDEsMi4xMi0uMDgsMi40LTUuMDcsNC45My04LjMxLDQuMzktNS42NCwxNC4zOS00LjU0LDE1LjMzLTQuNDIsMCwwLDMuMzEuNDEsMTYuMzcsNy43OWExMi42OCwxMi42OCwwLDAsMCwzLjYzLDIuNmM4LjM0LDMuOTIsMTQuODEtNC43OCwyNy4yOC00LDIuNDYuMTUuMjUuMzcsMTUuODUsNSw4LjI2LDIuNDgsMTEuOCwzLjMsMTQsNi43NSwyLDMuMTUuNTcsNC43NSwzLjEyLDcuOHM1LjU1LDMsNy41Myw2LjQ5YTEwLjcyLDEwLjcyLDAsMCwxLDEuMyw0Ljk0Yy4zMywzLjkyLDguNzMsMTAuODcsNDIuODcsMjQuMTYiIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0yNi44MyAtMjEuNjcpIi8+PHBhdGggY2xhc3M9ImNscy0yIiBkPSJNMzQ1LDIzMS42N2ExOS44MSwxOS44MSwwLDAsMCw1LjE1LTYuMzFjMi4zMy00LjU4LDEuMTItNi45NCwzLjk1LTEwLjY1YTEyLjM0LDEyLjM0LDAsMCwxLDYuMzEtNC41NGMxLjM3LS40LDMtMS4xNSw0LjE0LS4zOS44NS41Ni44NiwxLjQ5LDIsMmEyLjc1LDIuNzUsMCwwLDAsMi4xNywwYy44Ni0uNDQsMS4xNi0xLjM4LDEuNTgtMi43Ni42Mi0yLjA3LjI5LTIuNy43OS0zLjE2Ljk0LS44NywzLjUxLjEsNS41MiwxLjM4LDYuMzYsNC4wNiw3LjMzLDExLjM2LDkuMjcsMTEsMS43NC0uMjksMS42LTYuMzMsMy45NC02LjcxLDEuMDYtLjE3LDEuMjgsMSwyLjc2LDEuMTksMS45Mi4xOCwzLjgtMS44NSw3LjQ5LTUuOTJhMjYuNTIsMjYuNTIsMCwwLDEsMS43OC0yYzIuMDgtMS44Miw2LjEtMy45NSw4LjY4LTIuNTcsMS4yOC43LDEuMTMsMS42OCwzLjk0LDQuMTQuOTEuODEsMS4zOCwxLjIxLDEuNzcsMS4xOSwxLjYyLS4wOCwxLjgzLTMuNzUsNS4xMy04LjQ4LDEuNjQtMi4zNSwyLjc2LTMuMTMsMy45NC0zLjM1LDEuMzUtLjI2LDEuNjMuNDIsMywuMTlzMS42Ny0uOTMsNS4xMy00LjE0YzEuNjEtMS40OSwyLjI4LTIsMy4xNS0yLDEuMTUuMDksMS4zOSwxLjExLDIuMTcsMS4xOCwxLjg3LjE4LDIuNzQtNS4zNSw2LjEyLTkuODZzOS4zMS04LjE5LDEyLjYyLTYuN2MxLjMzLjYsMS44NCwxLjg4LDUuMTIsMy43NSwxLjI5LjczLDEuNzcuODUsMi4zNywxLjE4LDUuMTUsMi44NSw1LDEzLjQ3LDUuOTIsMTMuNDEuNDQsMCwuMTMtMi4zMiwxLjc3LTMuNzVhNCw0LDAsMCwxLDMtMWMyLjEyLjMzLDIuMTIsMi45Miw1LjEzLDQuOTMuNy40NywyLjg5LDEuOTMsNC4zMywxLjE4LDEuOTEtMSwxLjI1LTUuMiwyLjE3LTUuMzNzMS40NSwyLjg0LDIuMTcsMi43N2MxLjE5LS4xMy42OC04LjE5LDIuMTctOC40OHM0LjQyLDguMjEsNS40Niw3LjljLjY2LS4yLDAtMy43MS44NS00czEuOTUsMi42MSwzLDIuMzdjLjg1LS4yMS42My0yLjM4LDEuMzgtMi41NywxLjI4LS4zMywzLjEyLDUuNyw0LjkzLDUuNTIsMS4xNy0uMTEsMS4xLTMuNDgsMi0zLjU1czEuNzcsNC4yNSwzLjk0LDkuNDdjMS42Nyw0LDIuMzUsNi4wOCwzLjE2LDYuMzFDNTIwLDIwNS43NSw1MzAsMTkzLjE4LDUzNSwxNzljMS41NS00LjQ1LDIuMzQtOC4zNCw1LjcyLTkuODZzNy43Ni0uMiwxMi42MiwxLjE5YTI3LjExLDI3LjExLDAsMCwxLDUuNTIsMmMzLjA4LDEuNTksNC42NiwzLjM2LDUuNzIsMi43NnMuNi0zLjA2LDEuMTktMy4xNSwxLjQxLDMuMjEsMi4zNiwzLjE1Yy44LDAsLjkxLTIuNDQsMS43OC0yLjU2czEuNTEsMS45MSwyLjM2LDMuMzVjMi42OSw0LjU1LDYuNzEsMywxMSw3LjQ5LDIuNjYsMi43OSwxLjY2LDMuOTMsNC45Myw3LjMsMi42NiwyLjc1LDQuMzMsMyw1LjkxLDQuNTMsMi4yOCwyLjE3LDQuMTQsNi42NywxLjkzLDE3IiB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMjYuODMgLTIxLjY3KSIvPjwvc3ZnPg=="),
                 tags$a(href="https://rdds.it", "Laboratory of RNA and Disease Data Science, University of Trento", style="color: currentcolor;")))
    ), 
    align = "center")
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 10 * 1024^3)
  # PROTN: Visibility of the enrichment parameters based on the value of the Enrichment radiobutton
  output$input_enrichment <- renderUI({
    if (input$enrichR) {
      tagList(
        fluidRow(
          radioButtons("enrichR_universe", "Execute enrichment of the whole Universe", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
          uiOutput("help21")
        ),
        fluidRow(
          selectizeInput("DB_enrich", "DB to analyse",
                         choices = lapply(split(read_tsv("R/dbs_enrichR.txt", col_names = FALSE)$X1,read_tsv("R/dbs_enrichR.txt", col_names = FALSE)[,2]), as.list),
                         selected = NULL, multiple = TRUE
          ),
          uiOutput("help15")
        ),
        fluidRow(
          radioButtons("enrichR_DB", "Execute enrichment only on selected database", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
          uiOutput("help22")
        ),
        fluidRow(
          textInput("terms_enrich", "Terms to search"),
          uiOutput("help18")
        ),
        fluidRow(
          radioButtons("pval_fdr_enrich", "Use FDR instead of P.Value for enrichment", c("TRUE", "FALSE"), inline = TRUE, selected = TRUE),
          uiOutput("help23")
        ),
        fluidRow(
          textInput("pvalue_enrich", "P.Value thr for enrichment", value = "0.05"),
          uiOutput("help16")
        ),
        fluidRow(
          sliderInput("os_enrich", "Overlap size thr for enrichment", 1, 30, step = 1, value = 5),
          uiOutput("help17")
        )
      )
    }
  })
  
  #PROTN: Visibility of the extra parameters based on the radiobutton
  output$input_params <- renderUI({
    if (input$custom_param){
      tagList(
        fluidRow(
          tags$hr()
        ),
        fluidRow(
          column(
            width = 4,
            fluidRow(
              textInput("filt_absent_value", "Threshold of acceptable missing values for condition", value = "0"),
              uiOutput("help19")
            ),
            fluidRow(
              radioButtons("pval_fdr", "Use FDR instead of P.Value", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
              uiOutput("help20")
            ),
            fluidRow(
              textInput("FC_DEPs", "Log2 FC thr", value = "0.75"),
              uiOutput("help8")
            ),
            fluidRow(
              textInput("pvalue_DEPs", "P.Value thr", value = "0.05"),
              uiOutput("help9")
            ),
            fluidRow(
              textInput("signal_DEPs", "Signal log2 expr thr", value = "-inf"),
              uiOutput("help7")
            )
          ),
          column(
            width = 4,
            fluidRow(
              selectizeInput("taxonomy", "NCBI Taxonomy ID", 
                             choice = data.table::fread("R/NCBI_taxID/subset_tax.csv", select = "name"), 
                             selected = "Homo sapiens", multiple = F),
              # uiOut.put("help10")
            ),
            fluidRow(
              radioButtons("batch_corr", "Batch effect correction", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
              uiOutput("help10")
            ),
            fluidRow(
              textInput("prot_boxplot", "Control Boxplot proteins"),
              uiOutput("help11")
            ),
            fluidRow(
              radioButtons("STRING", "Execute PPI network STRINGdb", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
              uiOutput("help13")
            ),
            fluidRow(
              radioButtons("enrichR", "Execute enrichment", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
              uiOutput("help14")
            )
          ),
          column(
            width = 4,
            uiOutput("input_enrichment")
          )
        )
      )
    }
  })
  
  #PHOSPROTN: Visibility of the enrichment parameters based on the value of the Enrichment radiobutton
  output$input_enrichment_phos <- renderUI({
    if (input$enrichR_phos) {
      tagList(
        fluidRow(
          radioButtons("enrichR_universe_phos", "Execute enrichment of the whole Universe", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
          uiOutput("help21_P")
        ),
        fluidRow(
          selectizeInput("DB_enrich_phos", "DB to analyse",
                         choices = lapply(split(read_tsv("R/dbs_enrichR.txt", col_names = FALSE)$X1,read_tsv("R/dbs_enrichR.txt", col_names = FALSE)[,2]), as.list),
                         selected = NULL, multiple = TRUE
          ),
          uiOutput("help18_P")
        ),
        fluidRow(
          radioButtons("enrichR_DB_phos", "Execute enrichment only on selected database", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
          uiOutput("help22_P")
        ),
        fluidRow(
          textInput("terms_enrich_phos", "Terms to search"),
          uiOutput("help15_P")
        ),
        fluidRow(
          radioButtons("pval_fdr_enrich_phos", "Use FDR instead of P.Value for enrichment", c("TRUE", "FALSE"), inline = TRUE, selected = TRUE),
          uiOutput("help23_P")
        ),
        fluidRow(
          textInput("pvalue_enrich_phos", "P.Value thr for enrichment", value = "0.05"),
          uiOutput("help16_P")
        ),
        fluidRow(
          sliderInput("os_enrich_phos", "Overlap size thr for enrichment", 1, 30, step = 1, value = 5),
          uiOutput("help17_P")
        )
      )
    }
  })
  
  #PHOSPROTN: Visibility of the extra parameters based on the radiobutton
  output$input_params_phos <- renderUI({
    if (input$custom_param_phos){
      tagList(
        fluidRow(
          tags$hr()
        ),
        fluidRow(
          column(
            width = 4,
            fluidRow(
              textInput("filt_absent_value_phos", "Threshold of acceptable missing values for condition", value = "0"),
              uiOutput("help19_P")
            ),
            fluidRow(
              radioButtons("pval_fdr_phos", "Use FDR instead of P.Value", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
              uiOutput("help20_P")
            ),
            fluidRow(
              textInput("phospho_phos", "Phosphorilation accuracy percentage thr (%)", value = "75"),
              uiOutput("help24_P")
            ),
            fluidRow(
              textInput("signal_DEPs_phos", "Signal log2 expr thr", value = "-inf"),
              uiOutput("help7_P")
            ),
            fluidRow(
              textInput("FC_DEPs_phos", "Log2 FC thr", value = "0.75"),
              uiOutput("help8_P")
            ),
            fluidRow(
              textInput("pvalue_DEPs_phos", "P.Value thr", value = "0.05"),
              uiOutput("help9_P")
            )
          ),
          column(
            width = 4,
            fluidRow(
              selectizeInput("taxonomy_phos", "NCBI Taxonomy ID", 
                             choice = data.table::fread("R/NCBI_taxID/subset_tax.csv", select = "name"), 
                             selected = "Homo sapiens", multiple = F),
              # uiOut.put("help10")
            ),
            fluidRow(
              radioButtons("batch_corr_phos", "Batch effect correction", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
              uiOutput("help10_P")
            ),
            fluidRow(
              textInput("prot_boxplot_phos", "Control Boxplot proteins"),
              uiOutput("help11_P")
            ),
            fluidRow(
              radioButtons("STRING_phos", "Execute PPI network STRINGdb", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
              uiOutput("help13_P")
            ),
            fluidRow(
              radioButtons("enrichR_phos", "Execute enrichment", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
              uiOutput("help14_P")
            ),
            fluidRow(
              radioButtons("kinaseTree_phos", "Draw the kinase trees", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
              uiOutput("help25_P")
            )
          ),
          column(
            width = 4,
            uiOutput("input_enrichment_phos")
          )
        )
      )
    }
  })
  
  #PHOSPROTN: Visibility of the proteomics files for PhosProTN
  output$input_filePROT_phos <- renderUI({
    if (input$sw_analyzer_phos == "ProteomeDiscoverer"){
      tagList(
        fluidRow(
          fileInput("input_file_prot", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
          uiOutput("help4_P")
        ),
        fluidRow(
          fileInput("pep_file_prot", "Select the PEP file of the PROTEOMICS..."),
          uiOutput("help5_P")
        ),
        fluidRow(
          fileInput("prot_file_prot", "Select the PROT file of the PROTEOMICS..."),
          uiOutput("help6_P")
        )
      )
    } else if(input$sw_analyzer_phos == "MaxQuant"){
      tagList(
        fluidRow(
          fileInput("input_file_prot", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
          uiOutput("help4_P")
        ),
        fluidRow(
          fileInput("pep_file_prot", "Select the EVIDENCE file of the PROTEOMICS..."),
          uiOutput("help26_P")
        )
      )
    } else{
      tagList(
        tags$p("BACK")
      )
    }
  })
  
  #PHOSPROTN: Visibility of the phospho-proteomics files for PhosProTN
  output$input_filePHOS_phos <- renderUI({
    if (input$sw_analyzer_phos == "ProteomeDiscoverer"){
      tagList(
        fluidRow(
          fileInput("input_file_phos", "Select the SAMPLE_ANNOTATION file of the PHOSPHOproteomics..."),
          uiOutput("help27_P")
        ),
        fluidRow(
          fileInput("pep_file_phos", "Select the PEP file of the PHOSPHOproteomics..."),
          uiOutput("help28_P")
        ),
        fluidRow(
          fileInput("prot_file_phos", "Select the PROT file of the PHOSPHOproteomics..."),
          uiOutput("help29_P")
        ),
        fluidRow(
          fileInput("psm_file_phos", "Select the PSM file of the PHOSPHOproteomics..."),
          uiOutput("help30_P")
        ),
        radioButtons("custom_param_phos", "Use custom parameter", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
      )
    } else if(input$sw_analyzer_phos == "MaxQuant"){
      tagList(
        fluidRow(
          fileInput("input_file_phos", "Select the SAMPLE_ANNOTATION file of the PHOSPHOproteomics..."),
          uiOutput("help27_P")
        ),
        fluidRow(
          fileInput("pep_file_phos", "Select the EVIDENCE file of the PHOSPHOproteomics..."),
          uiOutput("help31_P")
        ),
        radioButtons("custom_param_phos", "Use custom parameter", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
      )
    } else{
      tagList(
        tags$p("BACK")
      )
    }
  })
  
  #PROTN: Show modals with example files
  output$input_df <- DT::renderDataTable({
    input_file <- readxl::read_xlsx("Data/proteome/Input.xlsx")
    DT::datatable(input_file, escape = FALSE)
    
  })
  
  output$design_df <- DT::renderDataTable({
    design_file <- readxl::read_xlsx("Data/proteome/design.xlsx")
    DT::datatable(design_file, escape = FALSE)
  })
  
  output$pep_df <- DT::renderDataTable({
    pep_file <- data.table::fread("Data/proteome/peptides.txt", nrows = 100)
    DT::datatable(data.table::as.data.table(lapply(pep_file, function(x){substring(x, 1, 20)})), escape = FALSE, options=list(scrollX = T))
    
  })
  
  output$prot_df <- DT::renderDataTable({
    prot_file <- data.table::fread("Data/proteome/proteinGroups.txt", nrows = 100)
    DT::datatable(data.table::as.data.table(lapply(prot_file, function(x){substring(x, 1, 20)})), escape = FALSE, options=list(scrollX = T))
  })

  #PHOSPROTN: Show modals with example files phospho
  output$input_df_phos <- DT::renderDataTable({
    input_file_phos <- readxl::read_xlsx("Data/proteome/Input.xlsx")
    DT::datatable(input_file_phos, escape = FALSE)
    
  })
  
  output$design_df_phos <- DT::renderDataTable({
    design_file_phos <- readxl::read_xlsx("Data/proteome/design.xlsx")
    DT::datatable(design_file_phos, escape = FALSE)
  })
  
  output$pep_df_phos <- DT::renderDataTable({
    pep_file_phos <- data.table::fread("Data/proteome/peptides.txt", nrows = 100)
    DT::datatable(data.table::as.data.table(lapply(pep_file_phos, function(x){substring(x, 1, 20)})), escape = FALSE, options=list(scrollX = T))
  })
  
  output$prot_df_phos <- DT::renderDataTable({
    prot_file_phos <- data.table::fread("Data/proteome/proteinGroups.txt", nrows = 100)
    DT::datatable(data.table::as.data.table(lapply(prot_file_phos, function(x){substring(x, 1, 20)})), escape = FALSE, options=list(scrollX = T))
  })
  
  
  #Define HELP button in execution pages
  setHelp <- function(x) {
    showNotification(
      type = "message",
      duration = 7,
      switch(x,
             "help_btn1" = "Title of the experiment. It will be the title of the web page report.",
             "help_btn2" = "Description of the current experiment. It is the first paragraph of the report.",
             "help_btn3" = "Determine with software was use to identify peptides and proteins. Choice: Protein Discoverer; MaxQuant; TMT with Protein Discoverer.",
             "help_btn4" = "File with the information about the samples and the correlation between replicate ID and condition  (WARNING: Condition name MUST contain at least 1 character!). [FURTHER DETAILS IN THE INFO TAB]",
             "help_btn5" = "Raw file of peptides obtained from PD or MQ (file peptides.txt).",
             "help_btn6" = "Raw file of protein groups obtained from PD or MQ (file proteinGroups.txt).",
             "help_btn7" = "Signal log2 expression threshold for the differential analysis. DEFAULT: DEPs if Signal_log2_thr > -∞ (No Limit, represent by value \"inf\" in the cell)",
             "help_btn8" = "Fold Change threshold for the differential analysis. DEFAULT: DEPs if log2_FC_thr = 0.75 (Up-regulated > 0.75, Down-regulated < -0.75)",
             "help_btn9" = "P.value threshold for the differential analysis. DEFAULT: DEPs if P_value_thr < 0.05",
             "help_btn10" = "Execution of the batch effect correction performed by proBatch. If TRUE, column MS_batch required in Input file.",
             "help_btn11" = "List of proteins used as control of the intensities. For each protein a boxplot is generated comparing the mean of the intensities group by condition. Separator: ,",
             "help_btn12" = "Excel file containing the formulas of the contrast comparison you want to analyse. [FURTHER DETAILS IN THE INFO TAB]",
             "help_btn14" = "Execution of the enrichment step.",
             "help_btn13" = "Execution of the network analysis.",
             "help_btn15" = "Write the word that you want to search in the results of EnrichR (EX: MYC, C-MYC, Senescence,...). If empty, no plots are returned.",
             "help_btn16" = "P.value threshold for the enrichment analysis. DEFAULT: Term is significant if P_value_Enrich_thr < 0.05.",
             "help_btn17" = "Overlap size threshold. The overlap size is the number of DEPs discovered in the enriched terms. DEFAULT: Term is significant if min_overlap_gene > 5.",
             "help_btn18" = "Write the DBs that you want to see in your plots. If empty, no plots are returned.",
             "help_btn19" = "The number correspond to the upper limit of acceptable NA in the condition. If 0, a peptide require at least a condition with all the intensities.",
             "help_btn20" = "Choice on the usege of the FDR or the standard P.Value for the DEPs. DEFAULT: standard P.Value",
             "help_btn21" = "Execution of the global enrichment on the whole experimental universe. Enrichment of all proteins, not only DEPs.",
             "help_btn22" = "Execute the enrichment only on the datasets selected. By default the analysis is performed on all 98 datasets, provide more results but take longer time.",
             "help_btn23" = "Choice on the usege of the FDR or the standard P.Value for the Enrichment. DEFAULT: standard P.Value",
             
             "help_btn1_P" = "Title of the experiment. It will be the title of the web page report.",
             "help_btn2_P" = "Description of the current experiment. It is the first paragraph of the report.",
             "help_btn3_P" = "Determine with software was use to identify peptides and proteins. PD: Protein Discoverer; MQ: MaxQuant.",
             "help_btn4_P" = "File with the information about the samples and the correlation between replicate ID and condition of the proteomics (WARNING: Condition name MUST contain at least 1 character!). [FURTHER DETAILS IN THE INFO TAB]",
             "help_btn5_P" = "Raw file of peptides obtained from PD.",
             "help_btn26_P" = "Raw file of peptides obtained from MQ (file evidence.txt).",
             "help_btn6_P" = "Raw file of protein groups obtained from PD.",
             "help_btn7_P" = "Signal log2 expression threshold for the differential analysis. DEFAULT: DEPs if Signal_log2_thr > -∞ (No Limit, represent by value \"inf\" in the cell)",
             "help_btn8_P" = "Fold Change threshold for the differential analysis. DEFAULT: DEPs if log2_FC_thr = 0.75 (Up-regulated > 0.75, Down-regulated < -0.75)",
             "help_btn9_P" = "P.value threshold for the differential analysis. DEFAULT: DEPs if P_value_thr < 0.05",
             "help_btn10_P" = "Execution of the batch effect correction performed by proBatch. If TRUE, column MS_batch required in Input file.",
             "help_btn11_P" = "List of proteins used as control of the intensities. For each protein a boxplot is generated comparing the mean of the intensities group by condition. Separator: ,",
             "help_btn12_P" = "Excel file containing the formulas of the contrast comparison you want to analyse. [FURTHER DETAILS IN THE INFO TAB]",
             "help_btn14_P" = "Execution of the enrichment step.",
             "help_btn13_P" = "Execution of the network analysis.",
             "help_btn15_P" = "Write the word that you want to search in the results of EnrichR (EX: MYC, C-MYC, Senescence,...). If empty, no plots are returned.",
             "help_btn16_P" = "P.value threshold for the enrichment analysis. DEFAULT: Term is significant if P_value_Enrich_thr < 0.05.",
             "help_btn17_P" = "Overlap size threshold. The overlap size is the number of DEPs discovered in the enriched terms. DEFAULT: Term is significant if min_overlap_gene > 5.",
             "help_btn18_P" = "Write the DBs that you want to see in your plots. If empty, no plots are returned.",
             "help_btn19_P" = "The number correspond to the upper limit of acceptable NA in the condition. If 0, a peptide require at least a condition with all the intensities.",
             "help_btn20_P" = "Choice on the usege of the FDR or the standard P.Value for the DEPs. DEFAULT: standard P.Value",
             "help_btn21_P" = "Execution of the global enrichment on the whole experimental universe. Enrichment of all proteins, not only DEPs.",
             "help_btn22_P" = "Execute the enrichment only on the datasets selected. By default the analysis is performed on all 98 datasets, provide more results but take longer time.",
             "help_btn23_P" = "Choice on the usege of the FDR or the standard P.Value for the Enrichment. DEFAULT: standard P.Value",
             "help_btn24_P" = "Filter on the phosphorilation accuracy. Maintain peptides with % > value.",
             "help_btn25_P" = "Execution of the kinase activity tree analysis.",
             "help_btn27_P" = "File with the information about the samples and the correlation between replicate ID and condition of the phosphoproteomic (WARNING: Condition name MUST contain at least 1 character!). [FURTHER DETAILS IN THE INFO TAB]",
             "help_btn28_P" = "Raw file of peptides obtained from PD.",
             "help_btn31_P" = "Raw file of peptides obtained from MQ (file evidence.txt).",
             "help_btn29_P" = "Raw file of protein groups obtained from PD.",
             "help_btn30_P" = "Raw file of the PSM obtained from PD.",
      )
    )
  }
  # SHOW HELP NOTIFICATION PROTN
  {
    output$help1 <- renderUI(actionButton("help_btn1", "", icon = icon("circle-question")))
    output$help2 <- renderUI(actionButton("help_btn2", "", icon = icon("circle-question")))
    output$help3 <- renderUI(actionButton("help_btn3", "", icon = icon("circle-question")))
    output$help4 <- renderUI(actionButton("help_btn4", "", icon = icon("circle-question")))
    output$help5 <- renderUI(actionButton("help_btn5", "", icon = icon("circle-question")))
    output$help6 <- renderUI(actionButton("help_btn6", "", icon = icon("circle-question")))
    output$help7 <- renderUI(actionButton("help_btn7", "", icon = icon("circle-question")))
    output$help8 <- renderUI(actionButton("help_btn8", "", icon = icon("circle-question")))
    output$help9 <- renderUI(actionButton("help_btn9", "", icon = icon("circle-question")))
    output$help10 <- renderUI(actionButton("help_btn10", "", icon = icon("circle-question")))
    output$help11 <- renderUI(actionButton("help_btn11", "", icon = icon("circle-question")))
    output$help12 <- renderUI(actionButton("help_btn12", "", icon = icon("circle-question")))
    output$help13 <- renderUI(actionButton("help_btn13", "", icon = icon("circle-question")))
    output$help14 <- renderUI(actionButton("help_btn14", "", icon = icon("circle-question")))
    output$help15 <- renderUI(actionButton("help_btn15", "", icon = icon("circle-question")))
    output$help16 <- renderUI(actionButton("help_btn16", "", icon = icon("circle-question")))
    output$help17 <- renderUI(actionButton("help_btn17", "", icon = icon("circle-question")))
    output$help18 <- renderUI(actionButton("help_btn18", "", icon = icon("circle-question")))
    output$help19 <- renderUI(actionButton("help_btn19", "", icon = icon("circle-question")))
    output$help20 <- renderUI(actionButton("help_btn20", "", icon = icon("circle-question")))
    output$help21 <- renderUI(actionButton("help_btn21", "", icon = icon("circle-question")))
    output$help22 <- renderUI(actionButton("help_btn22", "", icon = icon("circle-question")))
    output$help23 <- renderUI(actionButton("help_btn23", "", icon = icon("circle-question")))
    paste0("help_btn", 1:23) %>%
      map(~ observeEvent(input[[.x]], {
        setHelp(.x)
        runjs(paste0('setTimeout(function(){$("#', str_remove(.x, "_btn"), '").append($("#shiny-notification-panel"))},0);'))
      }))
  }
  
  # SHOW HELP NOTIFICATION PhosProTN
  {
    output$help1_P <- renderUI(actionButton("help_btn1_P", "", icon = icon("circle-question")))
    output$help2_P <- renderUI(actionButton("help_btn2_P", "", icon = icon("circle-question")))
    output$help3_P <- renderUI(actionButton("help_btn3_P", "", icon = icon("circle-question")))
    output$help4_P <- renderUI(actionButton("help_btn4_P", "", icon = icon("circle-question")))
    output$help5_P <- renderUI(actionButton("help_btn5_P", "", icon = icon("circle-question")))
    output$help6_P <- renderUI(actionButton("help_btn6_P", "", icon = icon("circle-question")))
    output$help7_P <- renderUI(actionButton("help_btn7_P", "", icon = icon("circle-question")))
    output$help8_P <- renderUI(actionButton("help_btn8_P", "", icon = icon("circle-question")))
    output$help9_P <- renderUI(actionButton("help_btn9_P", "", icon = icon("circle-question")))
    output$help10_P <- renderUI(actionButton("help_btn10_P", "", icon = icon("circle-question")))
    output$help11_P <- renderUI(actionButton("help_btn11_P", "", icon = icon("circle-question")))
    output$help12_P <- renderUI(actionButton("help_btn12_P", "", icon = icon("circle-question")))
    output$help13_P <- renderUI(actionButton("help_btn13_P", "", icon = icon("circle-question")))
    output$help14_P <- renderUI(actionButton("help_btn14_P", "", icon = icon("circle-question")))
    output$help15_P <- renderUI(actionButton("help_btn15_P", "", icon = icon("circle-question")))
    output$help16_P <- renderUI(actionButton("help_btn16_P", "", icon = icon("circle-question")))
    output$help17_P <- renderUI(actionButton("help_btn17_P", "", icon = icon("circle-question")))
    output$help18_P <- renderUI(actionButton("help_btn18_P", "", icon = icon("circle-question")))
    output$help19_P <- renderUI(actionButton("help_btn19_P", "", icon = icon("circle-question")))
    output$help20_P <- renderUI(actionButton("help_btn20_P", "", icon = icon("circle-question")))
    output$help21_P <- renderUI(actionButton("help_btn21_P", "", icon = icon("circle-question")))
    output$help22_P <- renderUI(actionButton("help_btn22_P", "", icon = icon("circle-question")))
    output$help23_P <- renderUI(actionButton("help_btn23_P", "", icon = icon("circle-question")))
    output$help24_P <- renderUI(actionButton("help_btn24_P", "", icon = icon("circle-question")))
    output$help25_P <- renderUI(actionButton("help_btn25_P", "", icon = icon("circle-question")))
    output$help26_P <- renderUI(actionButton("help_btn26_P", "", icon = icon("circle-question")))
    output$help27_P <- renderUI(actionButton("help_btn27_P", "", icon = icon("circle-question")))
    output$help28_P <- renderUI(actionButton("help_btn28_P", "", icon = icon("circle-question")))
    output$help29_P <- renderUI(actionButton("help_btn29_P", "", icon = icon("circle-question")))
    output$help30_P <- renderUI(actionButton("help_btn30_P", "", icon = icon("circle-question")))
    output$help31_P <- renderUI(actionButton("help_btn31_P", "", icon = icon("circle-question")))
    paste0("help_btn", 1:31, "_P") %>%
      map(~ observeEvent(input[[.x]], {
        setHelp(.x)
        runjs(paste0('setTimeout(function(){$("#', str_remove(.x, "_btn"), '").append($("#shiny-notification-panel"))},0);'))
      }))
  }
  
  
  #####################################################
  addResourcePath("basedir", tempdir())
  output$preview_results <- renderUI({
    tryCatch(
      {
        withProgress(message = "Rendering, please wait!", {
          shinyjs::hide("modal_preview_output_ProTN")
          message(session$token)
          #Disable click page
          shinyjs::disable("report")
          shinyjs::disable("case_study")
          js$pageDisable("none")
          
          #Creation directory for the results
          dirOutput_2 <- tempdir()
          currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
          dirOutput_1 <- paste("/", currentTime, "/", sep = "")
          dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
          dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
          message(dirOutput_Server)
          #Create subfolder for the files
          dir.create(file.path(dirOutput_Server, "figures"), showWarnings = FALSE)
          dir.create(file.path(dirOutput_Server, "data"), showWarnings = FALSE)
          dir.create(file.path(dirOutput_Server, "tables"), showWarnings = FALSE)
          dir.create(file.path(paste0(dirOutput_Server, "figures"),"Expression"), showWarnings = FALSE)
          dir.create(file.path(paste0(dirOutput_Server, "figures"),"PCA_MDS"), showWarnings = FALSE)
          # Set up parameters to pass to Rmd document
          params <- list(
            doc_title = input$title_exp,
            description = input$description_exp,
            readPD_files = if (input$sw_analyzer == "ProteomeDiscoverer") {
              TRUE
            } else {
              FALSE
            },
            readMQ_files = if (input$sw_analyzer == "MaxQuant") {
              TRUE
            } else {
              FALSE
            },
            readTMT_files = if (input$sw_analyzer == "TMT_PD") {
              TRUE
            } else {
              FALSE
            },
            file_input = input$input_file$datapath,
            file_prot = input$prot_file$datapath,
            file_pep = input$pep_file$datapath,
            taxonomy = input$taxonomy,
            filt_absent_value = if(is.null(input$filt_absent_value)){"0"}else{input$filt_absent_value},
            pval_fdr = if(is.null(input$pval_fdr)){FALSE}else{input$pval_fdr},
            signal_thr = if(is.null(input$signal_DEPs)){"inf"}else{input$signal_DEPs},
            fc_thr = if(is.null(input$FC_DEPs)){"0.75"}else{input$FC_DEPs},
            pval_thr = if(is.null(input$pvalue_DEPs)){"0.05"}else{input$pvalue_DEPs},
            batch_corr_exe = if(is.null(input$batch_corr)){FALSE}else{input$batch_corr},
            contr_design = input$design$datapath,
            prot_boxplot = if(is.null(input$prot_boxplot)){""}else{input$prot_boxplot},
            run_enrich = if(is.null(input$enrichR)){FALSE}else{input$enrichR},
            run_enrich_universe = if(!(is.null(input$enrichR) | is.null(input$enrichR_universe))){if(input$enrichR==T){input$enrichR_universe}else{FALSE}}else{FALSE},
            run_STRING = if(is.null(input$STRING)){FALSE}else{input$STRING},
            pval_fdr_enrich = if(is.null(input$pval_fdr_enrich)){TRUE}else{input$pval_fdr_enrich},
            pval_enrich_thr = if(is.null(input$pvalue_enrich)){"0.05"}else{input$pvalue_enrich},
            overlap_size_enrich_thr = if(is.null(input$os_enrich)){as.integer(5)}else{input$os_enrich},
            enrich_filter_term = input$terms_enrich,
            enrich_filter_DBs = input$DB_enrich,
            enrichR_DB = if(is.null(input$enrichR)){FALSE}else{input$enrichR_DB},
            dirOutput = dirOutput_Server
          )
          
          if(is.null(params$file_input) | is.null(params$file_prot) | is.null(params$file_pep) | is.null(params$contr_design)){
            stop("Error: No file provide")
          }
          #Render the notebook for the analysis
          rmarkdown::render("R/pipeline_elaborate_PD_files.Rmd",
                            output_file = "protn_report.html",
                            output_dir = dirOutput_Server,
                            params = params,
                            envir = new.env(parent = globalenv())
          )
          
          #Reactivate click
          shinyjs::enable("report")
          shinyjs::enable("case_study")
          shinyjs::show("modal_preview_output_ProTN")
          js$pageDisable("all")
          
          #Save folder for the download
          readr::write_csv(data.frame("session"=session$token,
                                      "outdir"=dirOutput_Server),
                           file = paste0(tempdir(),"/outdir_log_ProTN.log"), append = T)
          
          tags$iframe(src = paste0("basedir", dirOutput_1,"/protn_report.html"), height = "100%", width = "100%", scrolling = "yes")
        })
      },
      error = function(e) {
        #Create error report and reactivate the click in the page
        showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
        shinyjs::enable("report")
        shinyjs::enable("case_study")
        shinyjs::show("modal_preview_output_ProTN")
        shinyjs::hide("download_report")
        js$pageDisable("all")
        html_text<-str_replace(read_file("R/error.html"), 
                               pattern = "The page you’re looking for doesn’t exist.</p>", 
                               replacement = paste0("Description:", e, "</p>"))
        write_file(html_text, file = paste0(tempdir(), "/error.html"))
        tags$iframe(src = "basedir/error.html", height = "100%", width = "100%", scrolling = "yes")
      }
    )
  })
  
  # EXECUTE ProTN
  output$download_report <- downloadHandler(
    filename = "results_ProTN.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Prepraring files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            #Save folder for the download
            logs<-readr::read_csv(file = paste0(tempdir(),"/outdir_log_ProTN.log"),col_names = F)
            dirOutput_Server_2<-as.list(logs[which(logs$X1==session$token),"X2"])
            oldwd <- getwd()
            setwd(dirOutput_Server_2[[length(dirOutput_Server_2)]])
            files2zip <- list.files("./", recursive = TRUE)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  
  #######################################################
  output$preview_results_CS <- renderUI({
    tryCatch(
      {
        withProgress(message = "Rendering, please wait!", {
          shinyjs::hide("modal_preview_case_study_ProTN")
          message(session$token)
          #Disable click page
          shinyjs::disable("report")
          shinyjs::disable("case_study")
          js$pageDisable("none")
          
          #Creation directory for the results
          dirOutput_2 <- tempdir()
          currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
          dirOutput_1 <- paste("/", currentTime, "/", sep = "")
          dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
          dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
          message(dirOutput_Server)
          #Create subfolder for the files
          dir.create(file.path(dirOutput_Server, "figures"), showWarnings = FALSE)
          dir.create(file.path(dirOutput_Server, "data"), showWarnings = FALSE)
          dir.create(file.path(dirOutput_Server, "tables"), showWarnings = FALSE)
          dir.create(file.path(paste0(dirOutput_Server, "figures"),"Expression"), showWarnings = FALSE)	
          dir.create(file.path(paste0(dirOutput_Server, "figures"),"PCA_MDS"), showWarnings = FALSE)
          
          # Set up parameters to pass to Rmd document
          params <- list(
            doc_title = "Example case study",
            description = "Example case study: Phosphoproteomics reveals that Parkinson’s disease kinase LRRK2 regulates a subset of Rab GTPases. PRIDE: PXD003071. \n \n DOI: 10.7554/eLife.12813, PubMed: 26824392, Des: Steger M, Tonelli F, Ito G, Davies P, Trost M, Vetter M, Wachter S, Lorentzen E, Duddy G, Wilson S, Baptista MA, Fiske BK, Fell MJ, Morrow JA, Reith AD, Alessi DR, Mann M. Phosphoproteomics reveals that Parkinson's disease kinase LRRK2 regulates a subset of Rab GTPases. Elife. 2016 Jan 29;5. pii: e12813",
            readPD_files = FALSE,
            readMQ_files = TRUE,
            readTMT_files = FALSE,
            file_input = "../Data/proteome/Input.xlsx",
            file_prot = "../Data/proteome/proteinGroups.txt",
            file_pep = "../Data/proteome/peptides.txt",
            taxonomy = "Homo sapiens",
            filt_absent_value = "0",
            pval_fdr = FALSE,
            signal_thr = "inf",
            fc_thr = "0.75",
            pval_thr = "0.05",
            batch_corr_exe = FALSE,
            contr_design = "../Data/proteome/design.xlsx",
            prot_boxplot = "Park8, Rab7L1, Rab29, Park16, Rab8A, Rab8, Rab10, Rab12",
            run_enrich = TRUE,
            run_enrich_universe = FALSE,
            run_STRING = TRUE,
            pval_fdr_enrich = TRUE,
            pval_enrich_thr = "0.05",
            overlap_size_enrich_thr = as.integer(5),
            enrich_filter_term = NULL,
            enrich_filter_DBs = c(
              "GO_Molecular_Function_2023",
              "GO_Cellular_Component_2023",
              "GO_Biological_Process_2023",
              "KEGG_2021_Human",
              "BioPlanet_2019",
              "MGI_Mammalian_Phenotype_Level_4_2019",
              "MSigDB_Hallmark_2020",
              "Jensen_COMPARTMENTS",
              "DisGeNET"
            ),
            enrichR_DB = FALSE,
            dirOutput = dirOutput_Server
          )
          if(is.null(params$file_input) | is.null(params$file_prot) | is.null(params$file_pep) | is.null(params$contr_design)){
            stop("Error: No file provide")
          }
          #Render the notebook for the analysis
          rmarkdown::render("R/pipeline_elaborate_PD_files.Rmd",
                            output_file = "protn_report.html",
                            output_dir = dirOutput_Server,
                            params = params,
                            envir = new.env(parent = globalenv())
          )
          
          #Reactivate click
          shinyjs::enable("report")
          shinyjs::enable("case_study")
          shinyjs::show("modal_preview_case_study_ProTN")
          js$pageDisable("all")
          
          #Save folder for the download
          readr::write_csv(data.frame("session"=session$token,
                                      "outdir"=dirOutput_Server),
                           file = paste0(tempdir(),"/outdir_log_ProTN.log"), append = T)
          
          tags$iframe(src = paste0("basedir", dirOutput_1,"/protn_report.html"), height = "100%", width = "100%", scrolling = "yes")
        })
      },
      error = function(e) {
        #Create error report and reactivate the click in the page
        showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
        shinyjs::enable("report")
        shinyjs::enable("case_study")
        shinyjs::show("modal_preview_case_study_ProTN")
        shinyjs::hide("download_report_CS")
        js$pageDisable("all")
        html_text<-str_replace(read_file("R/error.html"), 
                               pattern = "The page you’re looking for doesn’t exist.</p>", 
                               replacement = paste0("Description:", e, "</p>"))
        write_file(html_text, file = paste0(tempdir(), "/error.html"))
        
        tags$iframe(src = "basedir/error.html", height = "100%", width = "100%", scrolling = "yes")
      }
    )
  })
  
  # EXECUTE ProTN
  output$download_report_CS <- downloadHandler(
    filename = "results_case_study.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Prepraring files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            #Save folder for the download
            logs<-readr::read_csv(file = paste0(tempdir(),"/outdir_log_ProTN.log"),col_names = F)
            dirOutput_Server_2<-as.list(logs[which(logs$X1==session$token),"X2"])
            oldwd <- getwd()
            setwd(dirOutput_Server_2[[length(dirOutput_Server_2)]])
            files2zip <- list.files("./", recursive = TRUE)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  # EXECUTE PhosProTN
  output$preview_results_phos <- renderUI({
    tryCatch(
      {
        withProgress(message = "Rendering, please wait!", {
          shinyjs::hide("modal_preview_output_PhosProTN")
          message(session$token)
          #Disable click page
          shinyjs::disable("report_phos")
          js$pageDisable("none")
          
          #Creation directory for the results
          dirOutput_2 <- tempdir()
          currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
          dirOutput_1 <- paste("/", currentTime, "/", sep = "")
          dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
          dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
          dirOutput_Server_2<-dirOutput_Server
          message(dirOutput_Server)
          #Create subfolder for the files
          dir.create(file.path(dirOutput_Server, "figures"), showWarnings = FALSE)
          dir.create(file.path(dirOutput_Server, "data"), showWarnings = FALSE)
          dir.create(file.path(dirOutput_Server, "tables"), showWarnings = FALSE)
          dir.create(file.path(paste0(dirOutput_Server, "figures"),"Expression"), showWarnings = FALSE)
          dir.create(file.path(paste0(dirOutput_Server, "figures"),"PCA_MDS"), showWarnings = FALSE)
          # Set up parameters to pass to Rmd document
          params <- list(
            doc_title = input$title_exp_phos,
            description = input$description_exp_phos,
            readPD_files = if (input$sw_analyzer_phos == "ProteomeDiscoverer") {
              TRUE
            } else {
              FALSE
            },
            readMQ_files = if (input$sw_analyzer_phos == "MaxQuant") {
              TRUE
            } else {
              FALSE
            },
            file_input_prot = input$input_file_prot$datapath,
            file_prot_prot = if(input$sw_analyzer_phos == "ProteomeDiscoverer"){input$prot_file_prot$datapath}else{NA},
            file_pep_prot = input$pep_file_prot$datapath,
            file_input_phos = input$input_file_phos$datapath,
            file_prot_phos = if(input$sw_analyzer_phos == "ProteomeDiscoverer"){input$prot_file_phos$datapath}else{NA},
            file_pep_phos = input$pep_file_phos$datapath,
            file_psm_phos = if(input$sw_analyzer_phos == "ProteomeDiscoverer"){input$psm_file_phos$datapath}else{NA},
            taxonomy = input$taxonomy_phos,
            filt_absent_value = if(is.null(input$filt_absent_value_phos)){"0"}else{input$filt_absent_value_phos},
            pval_fdr = if(is.null(input$pval_fdr_phos)){FALSE}else{input$pval_fdr_phos},
            phospho_thr = if(is.null(input$phospho_phos)){"75"}else{input$phospho_phos},
            signal_thr = if(is.null(input$signal_DEPs_phos)){"inf"}else{input$signal_DEPs_phos},
            fc_thr = if(is.null(input$FC_DEPs_phos)){"0.75"}else{input$FC_DEPs_phos},
            pval_thr = if(is.null(input$pvalue_DEPs_phos)){"0.05"}else{input$pvalue_DEPs_phos},
            batch_corr_exe = if(is.null(input$batch_corr_phos)){FALSE}else{input$batch_corr_phos},
            contr_design = input$design_phos$datapath,
            prot_boxplot = if(is.null(input$prot_boxplot_phos)){""}else{input$prot_boxplot_phos},
            run_enrich = if(is.null(input$enrichR_phos)){FALSE}else{input$enrichR_phos},
            run_enrich_universe = if(!(is.null(input$enrichR_phos) | is.null(input$enrichR_universe_phos))){if(input$enrichR_phos==T){input$enrichR_universe_phos}else{FALSE}}else{FALSE},
            run_STRING = if(is.null(input$STRING_phos)){FALSE}else{input$STRING_phos},
            run_kinaseTree = if(is.null(input$kinaseTree_phos)){FALSE}else{input$kinaseTree_phos},
            pval_fdr_enrich = if(is.null(input$pval_fdr_enrich_phos)){TRUE}else{input$pval_fdr_enrich_phos},
            pval_enrich_thr = if(is.null(input$pvalue_enrich_phos)){"0.05"}else{input$pvalue_enrich_phos},
            overlap_size_enrich_thr = if(is.null(input$os_enrich_phos)){as.integer(5)}else{input$os_enrich_phos},
            enrich_filter_term = input$terms_enrich_phos,
            enrich_filter_DBs = input$DB_enrich_phos,
            enrichR_DB = if(is.null(input$enrichR_phos)){FALSE}else{input$enrichR_DB_phos},
            dirOutput = dirOutput_Server
          )
          if(is.null(params$file_input) | is.null(params$file_prot) | is.null(params$file_pep) | is.null(params$contr_design)){
            stop("Error: No file provide")
          }
          #Render the notebook for the analysis
          rmarkdown::render("R/pipeline_elaborate_PD_file_PhosProTN.Rmd",
                            output_file = "phosprotn_report.html",
                            output_dir = dirOutput_Server,
                            params = params,
                            envir = new.env(parent = globalenv())
          )
          
          #Reactivate click
          shinyjs::enable("report_phos")
          js$pageDisable("all")
          shinyjs::show("modal_preview_output_PhosProTN")
          
          #Save folder for the download
          readr::write_csv(data.frame("session"=session$token,
                                      "outdir"=dirOutput_Server),
                           file = paste0(tempdir(),"/outdir_log_PhosProTN.log"), append = T)
          
          tags$iframe(src = paste0("basedir", dirOutput_1,"/phosprotn_report.html"), height = "100%", width = "100%", scrolling = "yes")
        })
      },
      error = function(e) {
        #Create error report and reactivate the click in the page
        showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
        shinyjs::enable("report_phos")
        shinyjs::show("modal_preview_output_PhosProTN")
        shinyjs::hide("download_report_phos")
        js$pageDisable("all")
        html_text<-str_replace(read_file("R/error.html"), 
                               pattern = "The page you’re looking for doesn’t exist.</p>", 
                               replacement = paste0("Description:", e, "</p>"))
        write_file(html_text, file = paste0(tempdir(), "/error.html"))
        
        
        tags$iframe(src = "basedir/error.html", height = "100%", width = "100%", scrolling = "yes")
      }
    )
  })
  
  # EXECUTE ProTN
  output$download_report_PhosProTN <- downloadHandler(
    filename = "results_PhosProTN.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Prepraring files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            #Save folder for the download
            logs<-readr::read_csv(file = paste0(tempdir(),"/outdir_log_PhosProTN.log"),col_names = F)
            dirOutput_Server_2<-as.list(logs[which(logs$X1==session$token),"X2"])
            oldwd <- getwd()
            setwd(dirOutput_Server_2[[length(dirOutput_Server_2)]])
            files2zip <- list.files("./", recursive = TRUE)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  # -- DELETE TEMP FILES WHEN SESSION ENDS -- #
  session$onSessionEnded(function() {
    if (dir.exists(tempdir())){unlink(list.files(tempdir(), full.names = T), recursive = T)}
  })
}

# Run the application
shinyApp(ui = ui, server = server, options = list(port = 8100))
