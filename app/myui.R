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
####### Variables ##########
############################

load("./InputFiles/01_ProteomeData.RData")
load("./InputFiles/02_TranscriptomeData.RData")

############################
######### Script ###########
############################

ui <-  navbarPage(title = HTML("<b>TtDDR</b>: <i><u>T</u>etrahymena <u>t</u>hermophila</i> <u>D</u>NA <u>D</u>amage <u>R</u>esponse data base"),
                  theme = shinytheme("lumen"),
                  windowTitle = "TtDDR",
                  tabPanel("Welcome",
                           fluidRow(column(12,
                                           h2("Abstract"),
                                           p("A tightly regulated DNA damage response is critical to the overall integrity of the genome. 
                                             Here, we combine transcriptomics and proteomics to study DNA damage response kinetics across 
                                             well-established treatments in the ciliate Tetrahymena thermophila. This extensive data set 
                                             of 6 conditions (HU, MMS, IR, HP, cisplatin and UV) and 7 time points (from 0 to 8 hours) 
                                             integrates over 250 transcriptome and proteome measurements. We observed upregulation of 
                                             known DNA repair proteins as well as a global dynamic response of not yet characterized 
                                             transcripts and proteins. Using artificial neural networks, we classify different expression
                                             trends in response to the damaging agents. These networks revealed both a core and specific
                                             global dynamic response to the different genotoxic stressors, highlighting unexpected pathway
                                             crosstalk. Ultimately, our study provides novel insights into the DNA damage response kinetics 
                                             in Tetrahymena."))),
                           fluidRow(column(12,
                                           h3("Experimental design & Overview"))),
                           fluidRow(column(2,
                                           br()),
                                    column(8,
                                           br(),
                                           img(src = "Image_01.png",  height = "100%", width = "100%"),
                                           br(),
                                           br(),
                                           br(),
                                           HTML("<b>Screen to explore kinetics of DNA damage response in Tetrahymena.</b>
                                                <b>A)</b> Schematic of the screen workflow. Cells were treated with a mutagenic agent, and samples 
                                                were harvested incrementally over eight hours. At each timepoint, samples were collected for 
                                                processing for both RNA-sequencing and quantitative mass spectrometry. 
                                                <b>B)</b> Venn diagram depicting the overlap between the identified transcripts and proteins. 
                                                <b>C)</b> Lollipop plot of enriched DNA damage repair factors.")),
                                    column(2,
                                           br())),
                           fluidRow(column(12,
                                           tags$hr(style="border-color: darkgrey;"))),
                           fluidRow(column(12,
                                           strong("A systems view on DNA damage response kinetics"),
                                           br(),
                                           br(),
                                           p("Emily Nischwitz",tags$sup("1,6"),",",
                                             "Vivien A.C. Schoonenberg",tags$sup("1,4,5,6"),",",
                                             "Rachel Mullner",tags$sup("1,2"),",",
                                             "Albert Fradera-Sola",tags$sup("1,2"),",",
                                             "Susanne Zimbelmann",tags$sup("1"),",",
                                             "Joshua J. Smith",tags$sup("3"),",",
                                             "Falk Butter",tags$sup("*,1"),
                                             style = "font-size:10pt;"),
                                           em(
                                             "1 Institute of Molecular Biology (IMB), 55128 Mainz, Germany,", tags$br(),
                                             "2 Institute of Molecular Virology and Cell Biology, Friedrich-Loeffler-Institut, 17493 Greifswald, Germany", tags$br(),
                                             "3 Department of Biomedical Science, Missouri State University, Springfield, MO, 65897, USA", tags$br(),
                                             "4 Current address:  Division of Hematology/Oncology, Boston Children’s Hospital, Harvard Medical School, Boston, MA, 02115, USA", tags$br(),
                                             "5 Current address: Molecular Pathology Unit & Center for Cancer Research, Massachusetts General Hospital and Harvard Medical School, Boston, MA, 02114, USA", tags$br(),
                                             "6 These authors contributed equally. Either author can be listed first.", tags$br(),
                                             "* Correspondence: falk.butter@fli.de (FB)", style = "font-size:9pt;"),
                                           br(),
                                           br(),
                                           HTML("<h style=font-size:10pt > App created by Albert Fradera-Sola in January 2025. <b>Last update on February 2025</b></h>"),
                                           p("Comments and bug reports to the following e-mail: ", mailtoR(email = "A.FraderaSola@imb-mainz.de",
                                                                                                           subject = "Comments and bugs: Tetrahymena DB shiny app",
                                                                                                           text = "A.FraderaSola@imb-mainz.de"), style = "font-size:9pt;")
                                           ))),
                  tabPanel(HTML("<h style=color:#8FBC8F >Proteome data</h>"),
                           fluidRow(column(12,
                                           HTML("<h3><b>Data selection</b></h3>"))),
                           fluidRow(column(12,
                                           p("Here you can select which data to work with. First (1) you can filter out your data per treatment and timepoint.
                                             Then (2), you can search for an ID of interest either using protein IDs or gene names. Once you find an interesting ID, you can click on the table
                                             to select. An overview of selected IDs (3), allows you to keep track of which IDs you are currently working with. If you delete an ID from the selection,
                                             you should click update to de-select it from the table. Finally (4), a table showing the protein values per treatment and timepoint is generated."))),
                           fluidRow(column(1,
                                           HTML("<h3><b>(1): Settings</b></h3>"),
                                           checkboxGroupInput(inputId = "proteome_ui_treatmentID",
                                                              label = h4("Select treatment:"),
                                                              choices = sort(unique(proteome_df_main$Treatment)),
                                                              selected = sort(unique(proteome_df_main$Treatment))),
                                           checkboxGroupInput(inputId = "proteome_ui_timepointID",
                                                              label = h4("Select timepoint:"),
                                                              choices = unique(proteome_df_main$Timepoint),
                                                              selected = unique(proteome_df_main$Timepoint))),
                                    column(3,
                                           HTML("<h3><b>(2): Select proteins:</b></h3>"),
                                           dataTableOutput("proteome_SelectTable"),
                                           HTML("<h3><b>(3): Selected proteins:</b></h3>"),
                                           selectizeInput('proteome_ui_selectedIDs',
                                                          NULL,
                                                          choices=NULL,
                                                          multiple=TRUE),
                                           actionButton("update_proteome_MainDataTable",
                                                        "Update selection")),
                                    column(8,
                                           HTML("<h3><b>(4): Protein values</b></h3>"),
                                           dataTableOutput("proteome_MainDataTable"))),
                           fluidRow(column(12,
                                           tags$hr(style="border-color: darkgrey;"))),
                  tabsetPanel(
                    tabPanel("Dynamicity Analysis",
                             fluidRow(column(6,
                                             HTML("<h4><b>Line plot</b></h4>")),
                                      column(6,
                                             HTML("<h4><b>Scatter plot</b></h4>"))),
                             fluidRow(column(6,
                                             HTML("<b>Treatment: </b>
                                                  <b style=color:#133E67 > CPT </b><b style=color:#EDC01C > HP </b>
                                                  <b style=color:#1EB54F > HU </b><b style=color:#CE534D > IR </b>
                                                  <b style=color:#8B4399 > MMS </b><b style=color:#1973BB > UV </b>")),
                                      column(6,
                                             bsCollapsePanel("Advanced settings",
                                                             column(6,
                                                                    sliderInput("proteome_Dynamicity", "Select a Dynamicity threshold:",
                                                                                min = -3.5, max = 0,
                                                                                value = -1.5, step = 0.5)),
                                                             column(6,
                                                                    sliderInput("proteome_LFQ", "Select an Intensity threshold:",
                                                                                min = 0, max = 5,
                                                                                value = 2, step = 0.5))))),
                             fluidRow(column(6,
                                             plotlyOutput("proteome_LinePlot",height = "1000px")),
                                      column(6,
                                             plotlyOutput("proteome_ScatterPlot",height = "1000px")))))),
                  tabPanel(HTML("<h style=color:#8B678B >Transcriptome data</h>"),
                           fluidRow(column(12,
                                           HTML("<h3><b>Data selection</b></h3>"))),
                           fluidRow(column(12,
                                           p("Here you can select which data to work with. First (1) you can filter out your data per treatment and timepoint.
                                             Then (2), you can search for an ID of interest either using protein IDs or gene names. Once you find an interesting ID, you can click on the table
                                             to select. An overview of selected IDs (3), allows you to keep track of which IDs you are currently working with. If you delete an ID from the selection,
                                             you should click update to de-select it from the table. Finally (4), a table showing the protein values per treatment and timepoint is generated."))),
                           fluidRow(column(1,
                                           HTML("<h3><b>(1): Settings</b></h3>"),
                                           checkboxGroupInput(inputId = "transcriptome_ui_treatmentID",
                                                              label = h4("Select treatment:"),
                                                              choices = sort(unique(transcriptome_df_main$Treatment)),
                                                              selected = sort(unique(transcriptome_df_main$Treatment))),
                                           checkboxGroupInput(inputId = "transcriptome_ui_timepointID",
                                                              label = h4("Select timepoint:"),
                                                              choices = unique(transcriptome_df_main$Timepoint),
                                                              selected = unique(transcriptome_df_main$Timepoint))),
                                    column(3,
                                           HTML("<h3><b>(2): Select proteins:</b></h3>"),
                                           dataTableOutput("transcriptome_SelectTable"),
                                           HTML("<h3><b>(3): Selected proteins:</b></h3>"),
                                           selectizeInput('transcriptome_ui_selectedIDs',
                                                          NULL,
                                                          choices=NULL,
                                                          multiple=TRUE),
                                           actionButton("update_transcriptome_MainDataTable",
                                                        "Update selection")),
                                    column(8,
                                           HTML("<h3><b>(4): Protein values</b></h3>"),
                                           dataTableOutput("transcriptome_MainDataTable"))),
                           fluidRow(column(12,
                                           tags$hr(style="border-color: darkgrey;"))),
                           tabsetPanel(
                             tabPanel("Dynamicity Analysis",
                                      fluidRow(column(6,
                                                      HTML("<h4><b>Line plot</b></h4>")),
                                               column(6,
                                                      HTML("<h4><b>Scatter plot</b></h4>"))),
                                      fluidRow(column(6,
                                                      HTML("<b>Treatment: </b>
                                                  <b style=color:#133E67 > CPT </b><b style=color:#EDC01C > HP </b>
                                                  <b style=color:#1EB54F > HU </b><b style=color:#CE534D > IR </b>
                                                  <b style=color:#8B4399 > MMS </b><b style=color:#1973BB > UV </b>")),
                                               column(6,
                                                      bsCollapsePanel("Advanced settings",
                                                                      column(6,
                                                                             sliderInput("transcriptome_Dynamicity", "Select a Dynamicity threshold:",
                                                                                         min = -3.5, max = 0,
                                                                                         value = -1.5, step = 0.5)),
                                                                      column(6,
                                                                             sliderInput("transcriptome_LFQ", "Select an Intensity threshold:",
                                                                                         min = 0, max = 5,
                                                                                         value = 2, step = 0.5))))),
                                      fluidRow(column(6,
                                                      plotlyOutput("transcriptome_LinePlot",height = "1000px")),
                                               column(6,
                                                      plotlyOutput("transcriptome_ScatterPlot",height = "1000px")))))))
