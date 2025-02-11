############################
####### Libraries ##########
############################

library(tidyverse)

############################
######## Script ############
############################

#### Proteome data

data_RAW <- read.csv(file = "Protein_data_calc.csv")

data <- data_RAW[data_RAW$score == "L2FC_LFQ.int",]

Encoding(data$Gene.Name) <- "UTF-8"

data$Gene.Name <- iconv(data$Gene.Name, "UTF-8", "UTF-8",sub='')

data <- data[,c(1:10,12)]

data <- data %>% pivot_longer(cols = starts_with("H"),
                              names_to = "Timepoint",
                              values_to = "Value")

data <- data[,c(1:3,5,6,4)]

data$gini <- log10(data$gini)

# data$Database <- sprintf('<a href="https://tet.ciliate.org/feature_details.php?feature_name=%s" target="_blank" class="btn btn-link" role="button">TGD</a>',data$Majority.protein.IDs)

data$Database <- sprintf('<a href="https://tet.ciliate.org/search.php?gene_name=%s" target="_blank" class="btn btn-link" role="button">TGD</a>',data$Majority.protein.IDs)

colnames(data)[3] <- "Treatment" 

proteome_df_main <- data

data <- data_RAW

data <- data_RAW[data_RAW$score == "L2FC_LFQ.int",]

Encoding(data$Gene.Name) <- "UTF-8"

data$Gene.Name <- iconv(data$Gene.Name, "UTF-8", "UTF-8",sub='')

data <- data[,grep(pattern = ".*protein.*|treatment|gini|big",x = colnames(data))]

data <- data %>% arrange(treatment,Majority.protein.IDs)

proteome_df_dynamicity <- data

save(proteome_df_main, proteome_df_dynamicity, file = "01_ProteomeData.RData")

### Transcriptome data

data_RAW <- read.csv(file = "RNA_data_calc.csv")

data <- data_RAW[data_RAW$score == "L2FC_CPM",]

Encoding(data$gene_name) <- "UTF-8"

data$gene_name <- iconv(data$gene_name, "UTF-8", "UTF-8",sub='')

data <- data[,c(1,2,4:11)]

data <- data %>% pivot_longer(cols = starts_with("H"),
                              names_to = "Timepoint",
                              values_to = "Value")

data <- data[,c(1:2,4,5,3)]

data$gini <- log10(data$gini)

# data$Database <- sprintf('<a href="https://tet.ciliate.org/feature_details.php?feature_name=%s" target="_blank" class="btn btn-link" role="button">TGD</a>',data$Majority.protein.IDs)

data$Database <- sprintf('<a href="https://tet.ciliate.org/search.php?gene_name=%s" target="_blank" class="btn btn-link" role="button">TGD</a>',data$gene_name)

colnames(data)[2] <- "Treatment" 

transcriptome_df_main <- data

data_RAW <- read.csv(file = "RNA_data_calc.csv")

data <- data_RAW[data_RAW$score == "L2FC_CPM",]

Encoding(data$gene_name) <- "UTF-8"

data$gene_name <- iconv(data$gene_name, "UTF-8", "UTF-8",sub='')

data <- data[,grep(pattern = ".*gene.*|treatment|gini|big",x = colnames(data))]

data <- data %>% arrange(treatment,gene_name)

transcriptome_df_dynamicity <- data

save(transcriptome_df_main, transcriptome_df_dynamicity, file = "02_TranscriptomeData.RData")
