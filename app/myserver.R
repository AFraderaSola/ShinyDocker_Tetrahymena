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
####### Functions ##########
############################

############################
######### Script ###########
############################

server <- function(input, output, session) {
  
  #### Data functions ####
  
  # Modify the different data sets in function of
  # user's input
  
  ##### Proteome ####
  
  ## Main data table
  
  proteome_SelectTable1 <- reactive({
    my_select_table_data <- proteome_df_main
    my_select_table_data <- my_select_table_data %>% select(Majority.protein.IDs, Gene.Name,Database) %>% distinct()
    colnames(my_select_table_data)[1] <- "Protein ID"
    colnames(my_select_table_data)[2] <- "Gene Name"
    my_select_table_data
  })
  
  ## Main data table
  
  proteome_MainDataTable1 <- reactive({
    my_main_table_data <- proteome_df_main
    my_main_table_data <- my_main_table_data %>% select(!Database)
    my_main_table_data <- my_main_table_data %>% arrange(Treatment, Majority.protein.IDs, Timepoint)
    if(length(input$proteome_ui_treatmentID) > 0) {
      my_main_table_data <- my_main_table_data %>% filter(Treatment %in% input$proteome_ui_treatmentID)
    }
    if(length(input$proteome_ui_timepointID) > 0) {
      my_main_table_data <- my_main_table_data %>% filter(Timepoint %in% input$proteome_ui_timepointID)
    }
    my_main_table_data
  })
  
  proteome_MainDataTable2 <- reactive({
    my_main_table_data <- proteome_MainDataTable1()
    colnames(my_main_table_data)[1] <- "Protein ID"
    colnames(my_main_table_data)[2] <- "Gene Name"
    colnames(my_main_table_data)[5] <- "FoldChange(log2[LFQ])"
    colnames(my_main_table_data)[6] <- "log10(Gini score)"
    my_main_table_data
  })
  
  ## Dynamic data table
  
  proteome_DynamicTable1 <- reactive({
    my_dynamic_table_data <- proteome_df_dynamicity
    if(length(input$proteome_ui_treatmentID) > 0) {
      my_dynamic_table_data <- my_dynamic_table_data %>% filter(treatment %in% input$proteome_ui_treatmentID)
    }
    my_dynamic_table_data
  })
  
  ##### Transcriptome ####
  
  ## Main data table
  
  transcriptome_SelectTable1 <- reactive({
    my_select_table_data <- transcriptome_df_main
    my_select_table_data <- my_select_table_data %>% select(gene_name,Database) %>% distinct()
    colnames(my_select_table_data)[1] <- "Gene Name"
    my_select_table_data
  })
  
  ## Main data table
  
  transcriptome_MainDataTable1 <- reactive({
    my_main_table_data <- transcriptome_df_main
    my_main_table_data <- my_main_table_data %>% select(!Database)
    my_main_table_data <- my_main_table_data %>% arrange(Treatment, gene_name, Timepoint)
    if(length(input$transcriptome_ui_treatmentID) > 0) {
      my_main_table_data <- my_main_table_data %>% filter(Treatment %in% input$transcriptome_ui_treatmentID)
    }
    if(length(input$transcriptome_ui_timepointID) > 0) {
      my_main_table_data <- my_main_table_data %>% filter(Timepoint %in% input$transcriptome_ui_timepointID)
    }
    my_main_table_data
  })
  
  transcriptome_MainDataTable2 <- reactive({
    my_main_table_data <- transcriptome_MainDataTable1()
    colnames(my_main_table_data)[1] <- "Gene Name"
    colnames(my_main_table_data)[4] <- "FoldChange(log2[CPM])"
    colnames(my_main_table_data)[5] <- "log10(Gini score)"
    my_main_table_data
  })
  
  ## Dynamic data table
  
  transcriptome_DynamicTable1 <- reactive({
    my_dynamic_table_data <- transcriptome_df_dynamicity
    if(length(input$transcriptome_ui_treatmentID) > 0) {
      my_dynamic_table_data <- my_dynamic_table_data %>% filter(treatment %in% input$transcriptome_ui_treatmentID)
    }
    my_dynamic_table_data
  })
  
  #### Output functions ####
  
  #### Proteome #### 
  
  #### Select table output  
  
  output$proteome_SelectTable <- renderDataTable({
    proteome_SelectTable1() %>%
      datatable(list(pageLength = 2, 
                     paging = T,    ## paginate the output
                     scrollX = T,   ## enable scrolling on X axis
                     scrollY = T,   ## enable scrolling on Y axis
                     autoWidth = T),
                rownames = F,
                escape = F)
  })
  
  proteome_SelectTable_proxy <- dataTableProxy('proteome_SelectTable')
  
  #### Main table output  
  
  output$proteome_MainDataTable <- renderDataTable({
    
    UI_ID <- input$proteome_ui_selectedIDs
    
    highlight <- unique(UI_ID)
    
    if (length(highlight) > 0) {
      
      my_main_data_table <- proteome_MainDataTable2()
      my_main_data_table <- my_main_data_table[my_main_data_table$`Protein ID` %in% highlight,]
      my_main_data_table <- my_main_data_table %>%
        datatable(options = list(pageLength = 7, 
                                 paging = T,    ## paginate the output
                                 scrollX = T,   ## enable scrolling on X axis
                                 scrollY = T,   ## enable scrolling on Y axis
                                 autoWidth = T,
                                 dom = 'Bfrtip',
                                 buttons = c('csv', 'excel')),
                  rownames = F,
                  escape = F,
                  extensions = 'Buttons',
                  selection = 'single', ## enable selection of a single row
                  filter = 'bottom') %>%
        formatRound("FoldChange(log2[LFQ])", digits = 2) %>%
        formatRound("log10(Gini score)", digits = 2)
      
      my_main_data_table
      
    }else{
      
      my_main_data_table <- proteome_MainDataTable2()
      my_main_data_table <- my_main_data_table %>%
        datatable(options = list(pageLength = 7, 
                                 paging = T,    ## paginate the output
                                 scrollX = T,   ## enable scrolling on X axis
                                 scrollY = T,   ## enable scrolling on Y axis
                                 autoWidth = T,
                                 dom = 'Bfrtip',
                                 buttons = c('csv', 'excel')),
                  rownames = F,
                  escape = F,
                  extensions = 'Buttons',
                  selection = 'single', ## enable selection of a single row
                  filter = 'bottom') %>%
        formatRound("FoldChange(log2[LFQ])", digits = 2) %>%
        formatRound("log10(Gini score)", digits = 2)
      
      my_main_data_table
      
    }
  })
  
  ##### Line Plot 
  
  output$proteome_LinePlot <- renderPlotly({
    
    if(length(input$proteome_ui_treatmentID) > 0){
      
      if (length(input$proteome_ui_timepointID) > 1) {
        
        # Proteins to highlight:
        
        UI_ID <- input$proteome_ui_selectedIDs
        
        highlight <- unique(UI_ID)
        
        if (length(highlight) > 0) {
          
          df <- proteome_MainDataTable1()
          
          # Retrieve the mean columns & tidy
          
          df <- df[,!grepl(pattern = "Gene.*|Database",x = colnames(df))]
          
          plot_df <- df
          
          plot_df <- plot_df[plot_df$Majority.protein.IDs %in% highlight,]
          
          plot_df$Treatment <- factor(x = plot_df$Treatment, levels = unique(plot_df$Treatment))
          
          plot_df$Timepoint <- factor(x = plot_df$Timepoint, levels = unique(plot_df$Timepoint))
          
          plot_df$Group <- interaction(plot_df$Majority.protein.IDs, plot_df$Treatment)
          
          plot_df$Hoover <- paste0(sprintf("Value: %.2f<br />Timepoint: %s<br />Protein ID: %s",
                                           plot_df$Value,plot_df$Timepoint, plot_df$Majority.protein.IDs))
          
          # plot_df$Interaction <- interaction(interaction(plot_df$Majority.protein.IDs, plot_df$Treatment), plot_df$Timepoint)
          
          colors_lineplot <- c("#133E67", "#EDC01C", "#1EB54F", "#CE534D", "#8B4399", "#1973BB")
          
          names(colors_lineplot) <- sort(unique(unique(proteome_df_main$Treatment)))
          
          colors_lineplot <- colors_lineplot[names(colors_lineplot) %in% sort(unique(plot_df$Treatment))]
          
          proteome_lineplot <- ggplot(plot_df, aes(x=Timepoint, y=Value, 
                                                   shape = Majority.protein.IDs, 
                                                   group = Group, color = Treatment,
                                                   text = Hoover))+
            facet_wrap(~Majority.protein.IDs)+
            geom_line(linewidth = 1)+
            geom_point(size = 4)+
            # geom_point(data = subset(plot_df, highlight == TRUE), color = "red")+
            labs(shape = "Protein ID")+
            ylab("FoldChange(log2[LFQ])")+
            xlab("Hours after treatment") +
            scale_color_manual(labels = unique(plot_df$Treatment),
                               values = colors_lineplot)+
            # scale_shape_manual(labels  = unique(plot_df$gene_name),
            #                    values = shapes_lineplot)+
            # scale_y_continuous(limits = c(24,32.5),breaks = seq(25, 32.5, by = 2.5))+
            # guides(color=FALSE)+
            theme_minimal()+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
            theme(legend.title = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 14),
                  legend.key.size = unit(1, 'cm'),
                  legend.position = "top",
                  strip.text = element_text(size = 20, face = "italic"))
          
          proteome_lineplot <- ggplotly(proteome_lineplot, tooltip='text')
          
          proteome_lineplot <- proteome_lineplot %>% layout(showlegend = FALSE)
          
          print(proteome_lineplot)
          
        }else{
          
          proteome_lineplot <- ggplot()+
            # scale_x_discrete(name = hours)+
            scale_y_continuous(limits = c(-1,1),breaks = seq(-1, 1, by = 0.50))+
            xlab("Hours after treatment")+
            ylab("FoldChange(log2[LFQ])")+
            theme_minimal()+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
            theme(legend.title = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 14),
                  legend.key.size = unit(1, 'cm'),
                  legend.position = "none")
          
          print(proteome_lineplot)
          
        }
        
      }else{
        
        proteome_lineplot <- ggplot()+
          # scale_x_discrete(name = hours)+
          scale_y_continuous(limits = c(-1,1),breaks = seq(-1, 1, by = 0.50))+
          xlab("Hours after treatment")+
          ylab("FoldChange(log2[LFQ])")+
          theme_minimal()+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 14),
                legend.key.size = unit(1, 'cm'),
                legend.position = "none")
        
        print(proteome_lineplot)
        
        showNotification(paste("Dynamicity analysis needs at least 2 timepoints!"), duration = 10,type = "error")
        
      }
    }else{
      
      proteome_lineplot <- ggplot()+
        # scale_x_discrete(name = hours)+
        scale_y_continuous(limits = c(-1,1),breaks = seq(-1, 1, by = 0.50))+
        xlab("Hours after treatment")+
        ylab("FoldChange(log2[LFQ])")+
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      print(proteome_lineplot)
      
      showNotification(paste("Please select at least 1 treatment!"), duration = 10,type = "error")

    }
  
  })
  
  output$proteome_ScatterPlot <- renderPlotly({
    
    if(length(input$proteome_ui_treatmentID) > 0){
      
      if (length(input$proteome_ui_timepointID) > 1) {
        
        # Proteins to highlight:
        
        UI_ID <- input$proteome_ui_selectedIDs
        
        highlight <- unique(UI_ID)
        
        if (length(highlight) > 0) {
        
          df <- proteome_DynamicTable1()
          
          df$gini <- log10(df$gini)
          
          df$highlight <- df$Majority.protein.IDs %in% highlight
          
          plot_df <- df
          
          plot_df$Dynamicity <- paste0(sprintf(" %.2f<br />Intensity: %.2f<br />Protein ID: %s",
                                           plot_df$gini,plot_df$big, plot_df$Majority.protein.IDs))

          # gini threshold
          
          GINI_THRES <- input$proteome_Dynamicity
          
          # LFQ threshold
          
          LFQ_THRES <- input$proteome_LFQ
          
          proteome_scatter_plot <- ggplot(plot_df, aes(x=big, y=gini, label =Dynamicity)) +
            facet_wrap(~treatment,scales = "free")+
            geom_hline(yintercept = GINI_THRES, linetype = "longdash", color = "#b4b4b4", linewidth = 1)+
            geom_vline(xintercept = LFQ_THRES, linetype = "longdash", color = "#b4b4b4", linewidth = 1)+
            geom_vline(xintercept = -LFQ_THRES, linetype = "longdash", color = "#b4b4b4", linewidth = 1)+
            geom_pointdensity(data = subset(plot_df, highlight == FALSE),
                              aes(x = big, y = gini),
                              # adjust = .1,
                              alpha = .3,inherit.aes = T)+
            geom_point(data = subset(plot_df, highlight == TRUE),
                       aes(x = big, y = gini),
                       color ="red", size = 3)+
            scale_color_viridis(option = "D"
                                # ,limits = c(0,1250),breaks = seq(0, 1250, by = 250)
                                )+
            # scale_y_continuous(limits = c(0,-3.5),breaks = seq(0, -3.5, by = -1))+
            # scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
            xlab("max(log2[FoldChange{LFQ}])")+
            ylab("log10(dynamicity)")+
            # annotate('text', label=paste(sum(y>GINI_THRES & x<LFQ_THRES), "proteins"), x=min(plot_df$LFQ) + 0.5, y=max(plot_df$Dynamicity) +0.05, size = 7) +
            # annotate('text', label=paste(sum(y<GINI_THRES & x<LFQ_THRES), "proteins"), x=min(plot_df$LFQ) + 0.5, y=min(plot_df$Dynamicity), size = 7) +
            # annotate('text', label=paste(sum(y>GINI_THRES & x>LFQ_THRES), "proteins"), x=max(plot_df$LFQ) - 0.5, y=max(plot_df$Dynamicity) +0.05, size = 7) +
            # annotate('text', label=paste(sum(y<GINI_THRES & x>LFQ_THRES), "proteins"), x=max(plot_df$LFQ) - 0.5, y=min(plot_df$Dynamicity) +0.05, size = 7) +
            theme_minimal()+
            theme(axis.text=element_text(size=18),
                  axis.title=element_text(size=20,face="bold"),
                  strip.text = element_text(size = 20, face = "bold"))
            # theme(legend.title = element_text(size = 16, face = "bold"),
            #       legend.text = element_text(size = 14),
            #       legend.key.size = unit(1, 'cm'),
            #       legend.position = "none")
          
          proteome_scatter_plot <- ggplotly(proteome_scatter_plot, tooltip = "label")
          
          print(proteome_scatter_plot)
          
        }else{
          
          proteome_scatter_plot <- ggplot()+
            scale_y_continuous(limits = c(0,-3.5),breaks = seq(0, -3.5, by = -1))+
            scale_x_continuous(limits = c(-5,5),breaks = seq(-5, 5, by = 2.5))+
            xlab("max(log2[FoldChange{LFQ}])")+
            ylab("log10(dynamicity)")+
            theme_minimal()+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
            theme(legend.title = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 14),
                  legend.key.size = unit(1, 'cm'),
                  legend.position = "none")
          
          print(proteome_scatter_plot)
          
        }
        
      }else{
        
        proteome_scatter_plot <- ggplot()+
          scale_y_continuous(limits = c(0,-3.5),breaks = seq(0, -3.5, by = -1))+
          scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
          xlab("max(log2[FoldChange{LFQ}])")+
          ylab("log10(dynamicity)")+
          theme_minimal()+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 14),
                legend.key.size = unit(1, 'cm'),
                legend.position = "none")
        
        print(proteome_scatter_plot)
        
        showNotification(paste("Dynamicity analysis needs at least 2 timepoints!"), duration = 10,type = "error")
        
      }
    }else{
      
      proteome_scatter_plot <- ggplot()+
        scale_y_continuous(limits = c(0,-3.5),breaks = seq(0, -3.5, by = -1))+
        scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
        xlab("max(log2[FoldChange{LFQ}])")+
        ylab("log10(dynamicity)")+
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      print(proteome_scatter_plot)
      
      showNotification(paste("Please select at least 1 treatment!"), duration = 10,type = "error")
      
    }
    
  })
  
  #### Transcriptome #### 
  
  #### Select table output  
  
  output$transcriptome_SelectTable <- renderDataTable({
    transcriptome_SelectTable1() %>%
      datatable(list(pageLength = 2, 
                     paging = T,    ## paginate the output
                     scrollX = T,   ## enable scrolling on X axis
                     scrollY = T,   ## enable scrolling on Y axis
                     autoWidth = T),
                rownames = F,
                escape = F)
  })
  
  transcriptome_SelectTable_proxy <- dataTableProxy('transcriptome_SelectTable')
  
  #### Main table output  
  
  output$transcriptome_MainDataTable <- renderDataTable({
    
    UI_ID <- input$transcriptome_ui_selectedIDs
    
    highlight <- unique(UI_ID)
    
    if (length(highlight) > 0) {
      
      my_main_data_table <- transcriptome_MainDataTable2()
      my_main_data_table <- my_main_data_table[my_main_data_table$`Gene Name` %in% highlight,]
      my_main_data_table <- my_main_data_table %>%
        datatable(options = list(pageLength = 7, 
                                 paging = T,    ## paginate the output
                                 scrollX = T,   ## enable scrolling on X axis
                                 scrollY = T,   ## enable scrolling on Y axis
                                 autoWidth = T,
                                 dom = 'Bfrtip',
                                 buttons = c('csv', 'excel')),
                  rownames = F,
                  escape = F,
                  extensions = 'Buttons',
                  selection = 'single', ## enable selection of a single row
                  filter = 'bottom') %>%
        formatRound("FoldChange(log2[CPM])", digits = 2) %>%
        formatRound("log10(Gini score)", digits = 2)
      
      my_main_data_table
      
    }else{
      
      my_main_data_table <- transcriptome_MainDataTable2()
      my_main_data_table <- my_main_data_table %>%
        datatable(options = list(pageLength = 7, 
                                 paging = T,    ## paginate the output
                                 scrollX = T,   ## enable scrolling on X axis
                                 scrollY = T,   ## enable scrolling on Y axis
                                 autoWidth = T,
                                 dom = 'Bfrtip',
                                 buttons = c('csv', 'excel')),
                  rownames = F,
                  escape = F,
                  extensions = 'Buttons',
                  selection = 'single', ## enable selection of a single row
                  filter = 'bottom') %>%
        formatRound("FoldChange(log2[CPM])", digits = 2) %>%
        formatRound("log10(Gini score)", digits = 2)
      
      my_main_data_table
      
    }
  })
  
  ##### Line Plot 
  
  output$transcriptome_LinePlot <- renderPlotly({
    
    if(length(input$transcriptome_ui_treatmentID) > 0){
      
      if (length(input$transcriptome_ui_timepointID) > 1) {
        
        # Proteins to highlight:
        
        UI_ID <- input$transcriptome_ui_selectedIDs
        
        highlight <- unique(UI_ID)
        
        if (length(highlight) > 0) {
          
          df <- transcriptome_MainDataTable1()
          
          # Retrieve the mean columns & tidy
          
          df <- df[,!grepl(pattern = "Database",x = colnames(df))]
          
          plot_df <- df
          
          plot_df <- plot_df[plot_df$gene_name %in% highlight,]
          
          plot_df$Treatment <- factor(x = plot_df$Treatment, levels = unique(plot_df$Treatment))
          
          plot_df$Timepoint <- factor(x = plot_df$Timepoint, levels = unique(plot_df$Timepoint))
          
          plot_df$Group <- interaction(plot_df$gene_name, plot_df$Treatment)
          
          plot_df$Hoover <- paste0(sprintf("Value: %.2f<br />Timepoint: %s<br />Protein ID: %s",
                                           plot_df$Value,plot_df$Timepoint, plot_df$gene_name))
          
          # plot_df$Interaction <- interaction(interaction(plot_df$gene_name, plot_df$Treatment), plot_df$Timepoint)
          
          colors_lineplot <- c("#133E67", "#EDC01C", "#1EB54F", "#CE534D", "#8B4399", "#1973BB")
          
          names(colors_lineplot) <- sort(unique(unique(transcriptome_df_main$Treatment)))
          
          colors_lineplot <- colors_lineplot[names(colors_lineplot) %in% sort(unique(plot_df$Treatment))]
          
          transcriptome_lineplot <- ggplot(plot_df, aes(x=Timepoint, y=Value, 
                                                   shape = gene_name, 
                                                   group = Group, color = Treatment,
                                                   text = Hoover))+
            facet_wrap(~gene_name)+
            geom_line(linewidth = 1)+
            geom_point(size = 4)+
            # geom_point(data = subset(plot_df, highlight == TRUE), color = "red")+
            labs(shape = "Protein ID")+
            ylab("FoldChange(log2[LFQ])")+
            xlab("Hours after treatment") +
            scale_color_manual(labels = unique(plot_df$Treatment),
                               values = colors_lineplot)+
            # scale_shape_manual(labels  = unique(plot_df$gene_name),
            #                    values = shapes_lineplot)+
            # scale_y_continuous(limits = c(24,32.5),breaks = seq(25, 32.5, by = 2.5))+
            # guides(color=FALSE)+
            theme_minimal()+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
            theme(legend.title = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 14),
                  legend.key.size = unit(1, 'cm'),
                  legend.position = "top",
                  strip.text = element_text(size = 20, face = "italic"))
          
          transcriptome_lineplot <- ggplotly(transcriptome_lineplot, tooltip='text')
          
          transcriptome_lineplot <- transcriptome_lineplot %>% layout(showlegend = FALSE)
          
          print(transcriptome_lineplot)
          
        }else{
          
          transcriptome_lineplot <- ggplot()+
            # scale_x_discrete(name = hours)+
            scale_y_continuous(limits = c(-1,1),breaks = seq(-1, 1, by = 0.50))+
            xlab("Hours after treatment")+
            ylab("FoldChange(log2[CPM])")+
            theme_minimal()+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
            theme(legend.title = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 14),
                  legend.key.size = unit(1, 'cm'),
                  legend.position = "none")
          
          print(transcriptome_lineplot)
          
        }
        
      }else{
        
        transcriptome_lineplot <- ggplot()+
          # scale_x_discrete(name = hours)+
          scale_y_continuous(limits = c(-1,1),breaks = seq(-1, 1, by = 0.50))+
          xlab("Hours after treatment")+
          ylab("FoldChange(log2[CPM])")+
          theme_minimal()+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 14),
                legend.key.size = unit(1, 'cm'),
                legend.position = "none")
        
        print(transcriptome_lineplot)
        
        showNotification(paste("Dynamicity analysis needs at least 2 timepoints!"), duration = 10,type = "error")
        
      }
    }else{
      
      transcriptome_lineplot <- ggplot()+
        # scale_x_discrete(name = hours)+
        scale_y_continuous(limits = c(-1,1),breaks = seq(-1, 1, by = 0.50))+
        xlab("Hours after treatment")+
        ylab("FoldChange(log2[CPM])")+
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      print(transcriptome_lineplot)
      
      showNotification(paste("Please select at least 1 treatment!"), duration = 10,type = "error")

    }
  
  })
  
  output$transcriptome_ScatterPlot <- renderPlotly({
    
    if(length(input$transcriptome_ui_treatmentID) > 0){
      
      if (length(input$transcriptome_ui_timepointID) > 1) {
        
        # Proteins to highlight:
        
        UI_ID <- input$transcriptome_ui_selectedIDs
        
        highlight <- unique(UI_ID)
        
        if (length(highlight) > 0) {
        
          df <- transcriptome_DynamicTable1()
          
          df$gini <- log10(df$gini)
          
          df$highlight <- df$gene_name %in% highlight
          
          plot_df <- df
          
          plot_df$Dynamicity <- paste0(sprintf(" %.2f<br />Intensity: %.2f<br />Protein ID: %s",
                                           plot_df$gini,plot_df$big, plot_df$gene_name))

          # gini threshold
          
          GINI_THRES <- input$transcriptome_Dynamicity
          
          # LFQ threshold
          
          LFQ_THRES <- input$transcriptome_LFQ
          
          transcriptome_scatter_plot <- ggplot(plot_df, aes(x=big, y=gini, label =Dynamicity)) +
            facet_wrap(~treatment,scales = "free")+
            geom_hline(yintercept = GINI_THRES, linetype = "longdash", color = "#b4b4b4", linewidth = 1)+
            geom_vline(xintercept = LFQ_THRES, linetype = "longdash", color = "#b4b4b4", linewidth = 1)+
            geom_vline(xintercept = -LFQ_THRES, linetype = "longdash", color = "#b4b4b4", linewidth = 1)+
            geom_pointdensity(data = subset(plot_df, highlight == FALSE),
                              aes(x = big, y = gini),
                              # adjust = .1,
                              alpha = .3,inherit.aes = T)+
            geom_point(data = subset(plot_df, highlight == TRUE),
                       aes(x = big, y = gini),
                       color ="red", size = 3)+
            scale_color_viridis(option = "D"
                                # ,limits = c(0,1250),breaks = seq(0, 1250, by = 250)
                                )+
            # scale_y_continuous(limits = c(0,-3.5),breaks = seq(0, -3.5, by = -1))+
            # scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
            xlab("max(log2[FoldChange{CPM}])")+
            ylab("log10(dynamicity)")+
            # annotate('text', label=paste(sum(y>GINI_THRES & x<LFQ_THRES), "proteins"), x=min(plot_df$LFQ) + 0.5, y=max(plot_df$Dynamicity) +0.05, size = 7) +
            # annotate('text', label=paste(sum(y<GINI_THRES & x<LFQ_THRES), "proteins"), x=min(plot_df$LFQ) + 0.5, y=min(plot_df$Dynamicity), size = 7) +
            # annotate('text', label=paste(sum(y>GINI_THRES & x>LFQ_THRES), "proteins"), x=max(plot_df$LFQ) - 0.5, y=max(plot_df$Dynamicity) +0.05, size = 7) +
            # annotate('text', label=paste(sum(y<GINI_THRES & x>LFQ_THRES), "proteins"), x=max(plot_df$LFQ) - 0.5, y=min(plot_df$Dynamicity) +0.05, size = 7) +
            theme_minimal()+
            theme(axis.text=element_text(size=18),
                  axis.title=element_text(size=20,face="bold"),
                  strip.text = element_text(size = 20, face = "bold"))
            # theme(legend.title = element_text(size = 16, face = "bold"),
            #       legend.text = element_text(size = 14),
            #       legend.key.size = unit(1, 'cm'),
            #       legend.position = "none")
          
          transcriptome_scatter_plot <- ggplotly(transcriptome_scatter_plot, tooltip = "label")
          
          print(transcriptome_scatter_plot)
          
        }else{
          
          transcriptome_scatter_plot <- ggplot()+
            scale_y_continuous(limits = c(0,-3.5),breaks = seq(0, -3.5, by = -1))+
            scale_x_continuous(limits = c(-5,5),breaks = seq(-5, 5, by = 2.5))+
            xlab("max(log2[FoldChange{CPM}])")+
            ylab("log10(dynamicity)")+
            theme_minimal()+
            theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
            theme(legend.title = element_text(size = 16, face = "bold"),
                  legend.text = element_text(size = 14),
                  legend.key.size = unit(1, 'cm'),
                  legend.position = "none")
          
          print(transcriptome_scatter_plot)
          
        }
        
      }else{
        
        transcriptome_scatter_plot <- ggplot()+
          scale_y_continuous(limits = c(0,-3.5),breaks = seq(0, -3.5, by = -1))+
          scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
          xlab("max(log2[FoldChange{CPM}])")+
          ylab("log10(dynamicity)")+
          theme_minimal()+
          theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
          theme(legend.title = element_text(size = 16, face = "bold"),
                legend.text = element_text(size = 14),
                legend.key.size = unit(1, 'cm'),
                legend.position = "none")
        
        print(transcriptome_scatter_plot)
        
        showNotification(paste("Dynamicity analysis needs at least 2 timepoints!"), duration = 10,type = "error")
        
      }
    }else{
      
      transcriptome_scatter_plot <- ggplot()+
        scale_y_continuous(limits = c(0,-3.5),breaks = seq(0, -3.5, by = -1))+
        scale_x_continuous(limits = c(24,38.5),breaks = seq(25, 38, by = 2.5))+
        xlab("max(log2[FoldChange{CPM}])")+
        ylab("log10(dynamicity)")+
        theme_minimal()+
        theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))+
        theme(legend.title = element_text(size = 16, face = "bold"),
              legend.text = element_text(size = 14),
              legend.key.size = unit(1, 'cm'),
              legend.position = "none")
      
      print(transcriptome_scatter_plot)
      
      showNotification(paste("Please select at least 1 treatment!"), duration = 10,type = "error")
      
    }
    
  })
  
  #### Observe functions ####
  
  ### Proteome
  
  observe({
    updateSelectizeInput(
      session,
      'proteome_ui_selectedIDs',
      choices=proteome_SelectTable1()[['Protein ID']][as.numeric(input$proteome_SelectTable_rows_selected)],
      selected=proteome_SelectTable1()[['Protein ID']][as.numeric(input$proteome_SelectTable_rows_selected)]
    )
  })

  observeEvent(input$update_proteome_SelectTable, {
    rows <-
      match(input$proteome_ui_selectedIDs,
            proteome_SelectTable1()[['Protein ID']])
    selectRows(proteome_SelectTable_proxy,
               as.numeric(rows))
  })
  
  ### Transcriptome
  
  observe({
    updateSelectizeInput(
      session,
      'transcriptome_ui_selectedIDs',
      choices=transcriptome_SelectTable1()[['Gene Name']][as.numeric(input$transcriptome_SelectTable_rows_selected)],
      selected=transcriptome_SelectTable1()[['Gene Name']][as.numeric(input$transcriptome_SelectTable_rows_selected)]
    )
  })
  
  observeEvent(input$update_transcriptome_SelectTable, {
    rows <-
      match(input$transcriptome_ui_selectedIDs,
            transcriptome_SelectTable1()[['Gene Name']])
    selectRows(transcriptome_SelectTable_proxy,
               as.numeric(rows))
  })
  
}
