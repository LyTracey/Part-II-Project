setwd("/Users/traceyly/Desktop/R Project/B12")


library(tidyverse)
library(readxl)
library(tidyverse)
library(metaMA) #rowVars
library(zoo) #rollamean
library(pracma) #movavg
library(microbenchmark) #microbenchmark
library(slider) #slide_dbl
library(data.table) #frollmean
library(plotly) #plotly
library(ggpubr)
library(RColorBrewer)
library(reshape2) #dcast
library(readxl)




# Importing dataframes
data_frame_full <- data.frame(readr::read_tsv("GSE112394_fpkms.tsv"), stringsAsFactors = FALSE)

B12_data <- data.frame(readr::read_csv("B12-RNA-seq_experiment_fpkm_values_AH.csv"), stringsAsFactors = FALSE)
enzymes <- data.frame(readr::read_csv("2019-11-11_C1-cycle-enzymes_AH.csv"), stringsAsFactors = FALSE)

#Rename columns

colnames(B12_data)[colnames(B12_data) == 'geneID'] <- 'gene_id'

colnames(enzymes)[colnames(enzymes) == 'GeneID'] <- 'gene_id'

#Separate chromosome column into chromosome name, and start and end positions

data_frame <- data_frame_full %>% 
  separate(locus, into = c("chromosome", "start", "end"), sep = ":|-", convert = TRUE)

data_frame <- data_frame %>% 
  separate(gene_id, into = c("gene_id", "annotation"), sep = -5, convert = TRUE)

B12_data <- B12_data %>% 
  separate(gene_id, into = c("gene_id", "annotation"), sep = -5, convert = TRUE)

data_frame <- data_frame %>% filter(str_detect(chromosome, "chromosome")) %>% arrange(chromosome)



#Join data frames

genome <- left_join(data_frame, B12_data, by = "gene_id")

enzymes <- left_join(enzymes, genome, by = "gene_id")

enzymes <- subset(enzymes, select = -annotation.y)

colnames(enzymes)[colnames(enzymes) == 'annotation.x'] <- 'annotation'

enzymes <- as.data.frame(enzymes)


#String variables
time_frame <- c("dark.11h",
                "dark.9h",
                "dark.7h",
                "dark.5h",
                "dark.3h",
                "dark.1h",
                "dark.0.5h",
                "light.0h",
                "light.0.5h",
                "light.1h",
                "light.3h",
                "light.5h",
                "light.7h",
                "light.9h",
                "light.11h",
                "dark.13h")

ancestral_m <- c("A_mB12_1", 
                 "A_mB12_2",
                 "A_mB12_3",
                 "A_mB12_5")

ancestral_p <- c("A_pB12_1", 
                 "A_pB12_2", 
                 "A_pB12_3", 
                 "A_pB12_5")

met_m <- c("M_mB12_1",
           "M_mB12_2",
           "M_mB12_3",
           "M_mB12_5")

met_p <- c("M_pB12_1",
           "M_pB12_2",
           "M_pB12_3",
           "M_pB12_5")


# Calculations 1 function

calculations <- function(data) {
  
  data_1 <- data
  
  #Adding a mean position column for each entry
  #data_1 <- data %>% mutate(mean_position = (data$start + data$end)/2)
  
  #Creating variance and mean columns time_frame
  time <- data_1 %>% select(contains(time_frame))
  data_1$time_var <- rowVars(time, na.rm = TRUE)
  data_1$mean_FPKM <- rowMeans(time, na.rm = TRUE)
  data_1$sum_FPKM <- apply(time, 1, sum)
  data_1$sd <- sqrt(data_1$time_var)
  
  #Creating variance and mean columns ancestral
  a_m <- data_1 %>% select(contains(ancestral_m))
  a_p <- data_1 %>% select(contains(ancestral_p))
  m_m <- data_1 %>% select(contains(met_m))
  m_p <- data_1 %>% select(contains(met_p))
  
  data_1$A_mB12_mean <- rowMeans(a_m, na.rm = TRUE)
  data_1$A_mB12_var <- rowVars(a_m, na.rm = TRUE)
  data_1$A_mB12_sd <- sqrt(data_1$A_mB12_var)
  
  data_1$A_pB12_mean <- rowMeans(a_p, na.rm = TRUE)
  data_1$A_pB12_var <- rowVars(a_p, na.rm = TRUE)
  data_1$A_pB12_sd <- sqrt(data_1$A_pB12_var)
  
  data_1$M_mB12_mean <- rowMeans(m_m, na.rm = TRUE)
  data_1$M_mB12_var <- rowVars(m_m, na.rm = TRUE)
  data_1$M_mB12_sd <- sqrt(data_1$M_mB12_var)
  
  
  data_1$M_pB12_mean <- rowMeans(m_p, na.rm = TRUE)
  data_1$M_pB12_var <- rowVars(m_p, na.rm = TRUE)
  data_1$M_pB12_sd <- sqrt(data_1$M_pB12_var)
  
  #remove time data
  # data_1 <- data_1 %>% select(!contains(time_frame))
  
  
  return(data_1)
}

a <- enzymes %>% as.data.frame()

a <- calculations(a)



#Log calculations - calculations II

log_calculations <- function(data) {
  
  data_1 <- data
  
  data_1$log_var <- data$time_var %>% log10()
  data_1$log_FPKM <- data$mean_FPKM %>% log10()
  data_1$log_sd <- data$sd %>% log10()
  data_1$log_sum <- data$sum_FPKM %>% log10()
  
  
  data_1$log_FPKM[data_1$log_FPKM == -Inf] <- -NA
  data_1$log_var[data_1$log_var == -Inf] <- -NA
  data_1$log_sd[data_1$log_sd == -Inf] <- -NA
  data_1[data_1 == -Inf] <- NA
  
  
  return(data_1)
}

a <- log_calculations(a)

  
  