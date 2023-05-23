#Install packages for Shiny App Launcher
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("shiny", dependencies = T)
install.packages("tidyverse", dependencies = T)
install.packages("markdown", dependencies = T)
install.packages("knitr", dependencies = T)
install.packages("shinydashboard", dependencies = T)
install.packages("shinydashboardPlus", dependencies = T)
install.packages("shinymaterial", dependencies = T)
install.packages("shinyjs", dependencies = T)
install.packages("magrittr", dependencies = T)
install.packages("dplyr", dependencies = T)
install.packages("stringr", dependencies = T)
install.packages("shinyBS", dependencies = T)
install.packages("rmdformats", dependencies = T)
install.packages("parallel", dependencies = T)
install.packages("doParallel", dependencies = T)
install.packages("foreach", dependencies = T)
install.packages("extrafont", dependencies = T)


#Install packages for ProTN
install.packages("devtools", dependencies = T)
install.packages("enrichR", dependencies = T)
install.packages("data.table", dependencies = T)
install.packages("ggplot2", dependencies = T)
install.packages("ggrepel", dependencies = T)
install.packages("igraph", dependencies = T)
install.packages("ggraph", dependencies = T)
install.packages("graphlayouts", dependencies = T)
install.packages("RColorBrewer", dependencies = T)
install.packages("ggfortify", dependencies = T)
install.packages("ggbeeswarm", dependencies = T)
install.packages("ggthemes", dependencies = T)
install.packages("tidyr", dependencies = T)
install.packages("scales", dependencies = T)
install.packages("qpdf", dependencies = T)
install.packages("corrplot", dependencies = T)
install.packages("lazyeval", dependencies = T)
install.packages("lubridate", dependencies = T)
install.packages("pheatmap", dependencies = T)
install.packages("reshape2", dependencies = T)
install.packages("readr", dependencies = T)
install.packages("rlang", dependencies = T)
install.packages("tibble", dependencies = T)
install.packages("wesanderson", dependencies = T)
install.packages("svgPanZoom", dependencies = T)
install.packages("plotly", dependencies = T)
install.packages("pdftools", dependencies = T)
install.packages("magick", dependencies = T)
install.packages("WGCNA", dependencies = T)

BiocManager::install("biomaRt")
BiocManager::install("DEqMS")
BiocManager::install("STRINGdb")
BiocManager::install("GO.db")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("pvca")
BiocManager::install("sva")

devtools::install_github("symbioticMe/proBatch", dependencies = T)
devtools::install_github("PYangLab/PhosR", dependencies = T)

#Install additional packages for PhosProTN
install.packages("data.tree", dependencies = T)
install.packages("rsvg", dependencies = T)
install.packages("shinyWidgets", dependencies = T)
install.packages("colourpicker", dependencies = T)