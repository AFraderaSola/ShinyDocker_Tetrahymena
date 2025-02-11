############################
######## Libraries #########
############################

library(shiny)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(plotly)
library(ggrepel)
library(kohonen)
library(effectsize)
library(tibble)
library(tidyr)
library(rMQanalysis)
library(shinyBS)
library(DT)
library(mailtoR)
library(viridis)
library(pheatmap)
library(ggpointdensity)
library(shinythemes)
set.seed(666)

############################
######### Script ###########
############################

options(stringsAsFactors = FALSE)

source('myui.R', local = TRUE)
source('myserver.R')

shinyApp(
  ui = ui,
  server = server
)

