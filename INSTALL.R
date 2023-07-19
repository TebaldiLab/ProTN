################################################################################
# ProTN: an integrative pipeline for complete analysis of proteomics           # 
# data from mass spectrometry                                                  #
# Laboratory of RNA and Disease Data Science, University of Trento             #
# Developer: Gabriele Tomè                                                     #
# Issue at: https://github.com/TebaldiLab/ProTN/issues                         #
# PI: Dr. Toma Tebaldi, PhD                                                    #
################################################################################
#Install packages for Shiny App Launcher
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("shiny", quietly = TRUE))
  install.packages("shiny", dependencies = T)
if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse", dependencies = T)
if (!require("markdown", quietly = TRUE))
  install.packages("markdown", dependencies = T)
if (!require("knitr", quietly = TRUE))
  install.packages("knitr", dependencies = T)
if (!require("shinydashboard", quietly = TRUE))
  install.packages("shinydashboard", dependencies = T)
if (!require("shinydashboardPlus", quietly = TRUE))
  install.packages("shinydashboardPlus", dependencies = T)
if (!require("shinymaterial", quietly = TRUE))
  install.packages("shinymaterial", dependencies = T)
if (!require("shinyjs", quietly = TRUE))
  install.packages("shinyjs", dependencies = T)
if (!require("magrittr", quietly = TRUE))
  install.packages("magrittr", dependencies = T)
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr", dependencies = T)
if (!require("stringr", quietly = TRUE))
  install.packages("stringr", dependencies = T)
if (!require("shinyBS", quietly = TRUE))
  install.packages("shinyBS", dependencies = T)
if (!require("rmdformats", quietly = TRUE))
  install.packages("rmdformats", dependencies = T)
if (!require("parallel", quietly = TRUE))
  install.packages("parallel", dependencies = T)
if (!require("doParallel", quietly = TRUE))
  install.packages("doParallel", dependencies = T)
if (!require("foreach", quietly = TRUE))
  install.packages("foreach", dependencies = T)
if (!require("extrafont", quietly = TRUE))
  install.packages("extrafont", dependencies = T)

#Install packages for ProTN
if (!require("devtools", quietly = TRUE))
  install.packages("devtools", dependencies = T)
if (!require("enrichR", quietly = TRUE))
  install.packages("enrichR", dependencies = T)
if (!require("data.table", quietly = TRUE))
  install.packages("data.table", dependencies = T)
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2", dependencies = T)
if (!require("ggrepel", quietly = TRUE))
  install.packages("ggrepel", dependencies = T)
if (!require("igraph", quietly = TRUE))
  install.packages("igraph", dependencies = T)
if (!require("ggraph", quietly = TRUE))
  install.packages("ggraph", dependencies = T)
if (!require("graphlayouts", quietly = TRUE))
  install.packages("graphlayouts", dependencies = T)
if (!require("RColorBrewer", quietly = TRUE))
  install.packages("RColorBrewer", dependencies = T)
if (!require("ggfortify", quietly = TRUE))
  install.packages("ggfortify", dependencies = T)
if (!require("ggbeeswarm", quietly = TRUE))
  install.packages("ggbeeswarm", dependencies = T)
if (!require("ggthemes", quietly = TRUE))
  install.packages("ggthemes", dependencies = T)
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr", dependencies = T)
if (!require("scales", quietly = TRUE))
  install.packages("scales", dependencies = T)
if (!require("qpdf", quietly = TRUE))
  install.packages("qpdf", dependencies = T)
if (!require("corrplot", quietly = TRUE))
  install.packages("corrplot", dependencies = T)
if (!require("lazyeval", quietly = TRUE))
  install.packages("lazyeval", dependencies = T)
if (!require("lubridate", quietly = TRUE))
  install.packages("lubridate", dependencies = T)
if (!require("pheatmap", quietly = TRUE))
  install.packages("pheatmap", dependencies = T)
if (!require("reshape2", quietly = TRUE))
  install.packages("reshape2", dependencies = T)
if (!require("readr", quietly = TRUE))
  install.packages("readr", dependencies = T)
if (!require("rlang", quietly = TRUE))
  install.packages("rlang", dependencies = T)
if (!require("tibble", quietly = TRUE))
  install.packages("tibble", dependencies = T)
if (!require("wesanderson", quietly = TRUE))
  install.packages("wesanderson", dependencies = T)
if (!require("svgPanZoom", quietly = TRUE))
  install.packages("svgPanZoom", dependencies = T)
if (!require("plotly", quietly = TRUE))
  install.packages("plotly", dependencies = T)
if (!require("pdftools", quietly = TRUE))
  install.packages("pdftools", dependencies = T)
if (!require("magick", quietly = TRUE))
  install.packages("magick", dependencies = T)
if (!require("WGCNA", quietly = TRUE))
  install.packages("WGCNA", dependencies = T)

if (!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")
if (!require("DEqMS", quietly = TRUE))
  BiocManager::install("DEqMS")
if (!require("STRINGdb", quietly = TRUE))
  BiocManager::install("STRINGdb")
if (!require("GO.db", quietly = TRUE))
  BiocManager::install("GO.db")
if (!require("impute", quietly = TRUE))
  BiocManager::install("impute")
if (!require("preprocessCore", quietly = TRUE))
  BiocManager::install("preprocessCore")
if (!require("pvca", quietly = TRUE))
  BiocManager::install("pvca")
if (!require("sva", quietly = TRUE))
  BiocManager::install("sva")

if (!require("proBatch", quietly = TRUE))
  devtools::install_github("symbioticMe/proBatch", dependencies = T)
if (!require("PhosR", quietly = TRUE))
  devtools::install_github("PYangLab/PhosR", dependencies = T)

#Install additional packages for PhosProTN
if (!require("data.tree", quietly = TRUE))
  install.packages("data.tree", dependencies = T)
if (!require("rsvg", quietly = TRUE))
  install.packages("rsvg", dependencies = T)
if (!require("shinyWidgets", quietly = TRUE))
  install.packages("shinyWidgets", dependencies = T)
if (!require("colourpicker", quietly = TRUE))
  install.packages("colourpicker", dependencies = T, )

#Check installation package
list.of.packages <- c("BiocManager","shiny","tidyverse","markdown","knitr",
                      "shinydashboard","shinydashboardPlus","shinymaterial",
                      "shinyjs","magrittr","dplyr","stringr","shinyBS",
                      "rmdformats","parallel","doParallel","foreach","extrafont",
                      "devtools","enrichR","data.table","ggplot2","ggrepel",
                      "igraph","ggraph","graphlayouts","RColorBrewer","ggfortify",
                      "ggbeeswarm","ggthemes","tidyr","scales","qpdf",
                      "corrplot","lazyeval","lubridate","pheatmap","reshape2",
                      "readr","rlang","tibble","wesanderson","svgPanZoom",
                      "plotly","pdftools","magick","WGCNA","biomaRt",
                      "DEqMS","STRINGdb","GO.db","impute","preprocessCore",
                      "pvca","sva","proBatch","PhosR","data.tree",
                      "rsvg","shinyWidgets","colourpicker")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  message(paste0("ERROR: Some error occur durig the installation of the current package:\n\t",
                 paste0(new.packages, collapse = "\n\t"),
                 "\nPlease try to install them again."))}
