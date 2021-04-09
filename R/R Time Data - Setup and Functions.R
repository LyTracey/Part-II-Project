#Not presented as R markdown becaus processing takes too long


setwd("~/Desktop/R Project")

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
library(MESS) #auc 
library(writexl) #write_xlsx
library(nortest) #lillie test
library(stringi) #stri_unique



# Importing dataframe
data_frame_full <- data.frame(readr::read_tsv("GSE112394_fpkms.tsv"), stringsAsFactors = FALSE)

#Separate chromosome column into chromosome name, and start and end positions
data_frame <- data_frame_full %>% 
  separate(locus, into = c("chromosome", "start", "end"), sep = ":|-", convert = TRUE)

data_frame <- data_frame %>% filter(str_detect(chromosome, "chromosome")) %>% arrange(chromosome)



#Create chromosome subsets

name <- function(x) {
  paste("chromosome_", x, sep = "")
}

names <- sapply(1:17, name)

split <- split(data_frame, data_frame$chromosome)
split <- split[names]

lapply(seq_along(split), function(x) {
  assign(names[x], split[[x]], envir=.GlobalEnv)
}
)


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


# Calculations 1 function

calculations <- function(data) {
  
  data_1 <- data
  
  #Adding a mean position column for each entry
  #data_1 <- data %>% mutate(mean_position = (data$start + data$end)/2)
  
  #Creating variance and mean columns
  time <- data_1 %>% select(contains(time_frame))
  data_1$time_var <- rowVars(time, na.rm = TRUE)
  data_1$mean_FPKM <- rowMeans(time, na.rm = TRUE)
  data_1$sum_FPKM <- apply(time, 1, sum)
  data_1$sd <- sqrt(data_1$time_var)
  data_1$CV <- data_1$sd / data_1$mean_FPKM
  
  
  #remove time data
  data_1 <- data_1 %>% select(!contains(time_frame))
  return(data_1)
}



#chromosome_1_expanded <- calculations(chromosome_1)

#Expanding dataframe function
expand_dataframe <- function(dff, start_i, end_i) {
  
  #Add gene length column
  dff_expanded <- dff %>% mutate(gene_length = ((dff[,end_i] - dff[,start_i])  + 1))
  
  #Expand data with replicates within range 
  dff_expanded <- as.data.frame(lapply(dff_expanded, rep, dff_expanded$gene_length))
  
  #Create ranges
  bp <- as.vector(mapply(FUN = ":", dff[,start_i] , dff[,end_i])) %>% unlist()
  
  #Add ranges
  dff_expanded <- dff_expanded %>% mutate(bp_position = bp)
  
  #create vector with full range
  dff_whole <- data.frame(bp_position = min(dff[,start_i]):max(dff[,end_i]))
  
  #full_join of
  joined <- full_join(dff_whole, dff_expanded, by= "bp_position")
  
  #join
  joined <- arrange(joined, bp_position)
  
  return(joined)
}



#Log calculations - calculations II


log_calculations <- function(data) {
  
  data_1 <- data
  
  #Calculating the log10
  data_1$log_FPKM <- data_1$mean_FPKM %>% log10()
  data_1$log_var <-  data_1$time_var %>% log10()
  data_1$log_sum <- data_1$sum_FPKM %>% log10()
  data_1$log_sd <- data_1$sd %>% log10()
  
  #Assigning NA to infinite values
  data_1[data_1$log_FPKM == -Inf, "log_FPKM"] <- -NA
  data_1[data_1$log_var == -Inf] <- -NA
  data_1[data_1$log_sd == -Inf] <- -NA
  
  #Assigning the averaged mean FPKM for the whole data to the infinite values
  data_1$log_FPKM[is.na(data_1$log_FPKM)] <- mean(data_1$log_FPKM, na.rm = TRUE)
  data_1$log_FPKM[is.na(data_1$log_FPKM)] <- mean(data_1$log_FPKM, na.rm = TRUE)
  
  return(data_1)
}



#Rolling average function 1

RA1 <- function(data, name, window) {
  data_1 <- data[!duplicated(data$bp_position),]
  data_1[,name] <- frollmean(data_1$log_FPKM, n = window, align = "center", na.rm = TRUE)
  return(data_1)
}



#Rolling average function 2

RA2 <- function(data, name, window) {
  data_1 <- data[!duplicated(data$bp_position, fromLast = TRUE),]
  data_1[,name] <- frollmean(data_1$log_FPKM, n = window, align = "center", na.rm = TRUE)
  return(data_1)
}




#Sampling function
dff_sample <- function(dff, start, end, window ) {
  s <- c(seq.int(from = start, to = end, by = window))
  data <- dff[c(s),]
  return(data)
}



#Remove non-constitutively expressed genes from unexpanded original nuclear dataframe

data_frame_2 <- data_frame
data_frame_2[data_frame_2 == 0] <- NA
data_frame_2 <- data_frame_2[complete.cases(data_frame_2), ]

data_frame_2 <- calculations(data_frame_2)
data_frame_2 <- log_calculations(data_frame_2)


##APPLYING FUNCTIONS


#Applying calculations 1
chromosome_1_expanded <- calculations(chromosome_1)
chromosome_2_expanded <- calculations(chromosome_2)
chromosome_3_expanded <- calculations(chromosome_3)
chromosome_4_expanded <- calculations(chromosome_4)
chromosome_5_expanded <- calculations(chromosome_5)
chromosome_6_expanded <- calculations(chromosome_6)
chromosome_7_expanded <- calculations(chromosome_7)
chromosome_8_expanded <- calculations(chromosome_8)
chromosome_9_expanded <- calculations(chromosome_9)
chromosome_10_expanded <- calculations(chromosome_10)
chromosome_11_expanded <- calculations(chromosome_11)
chromosome_12_expanded <- calculations(chromosome_12)
chromosome_13_expanded <- calculations(chromosome_13)
chromosome_14_expanded <- calculations(chromosome_14)
chromosome_15_expanded <- calculations(chromosome_15)
chromosome_16_expanded <- calculations(chromosome_16)
chromosome_17_expanded <- calculations(chromosome_17)


#Applying expand
chromosome_1_expanded <- expand_dataframe(chromosome_1_expanded, 4, 5)
chromosome_2_expanded <- expand_dataframe(chromosome_2_expanded, 4, 5)
chromosome_3_expanded <- expand_dataframe(chromosome_3_expanded, 4, 5)
chromosome_4_expanded <- expand_dataframe(chromosome_4_expanded, 4, 5)
chromosome_5_expanded <- expand_dataframe(chromosome_5_expanded, 4, 5)
chromosome_6_expanded <- expand_dataframe(chromosome_6_expanded, 4, 5)
chromosome_7_expanded <- expand_dataframe(chromosome_7_expanded, 4, 5)
chromosome_8_expanded <- expand_dataframe(chromosome_8_expanded, 4, 5)
chromosome_9_expanded <- expand_dataframe(chromosome_9_expanded, 4, 5)
chromosome_10_expanded <- expand_dataframe(chromosome_10_expanded, 4, 5)
chromosome_11_expanded <- expand_dataframe(chromosome_11_expanded, 4, 5)
chromosome_12_expanded <- expand_dataframe(chromosome_12_expanded, 4, 5)
chromosome_13_expanded <- expand_dataframe(chromosome_13_expanded, 4, 5)
chromosome_14_expanded <- expand_dataframe(chromosome_14_expanded, 4, 5)
chromosome_15_expanded <- expand_dataframe(chromosome_15_expanded, 4, 5)
chromosome_16_expanded <- expand_dataframe(chromosome_16_expanded, 4, 5)
chromosome_17_expanded <- expand_dataframe(chromosome_17_expanded, 4, 5)


#Applying log calculations

chromosome_1_expanded <- log_calculations(chromosome_1_expanded)
chromosome_2_expanded <- log_calculations(chromosome_2_expanded)
chromosome_3_expanded <- log_calculations(chromosome_3_expanded)
chromosome_4_expanded <- log_calculations(chromosome_4_expanded)
chromosome_5_expanded <- log_calculations(chromosome_5_expanded)
chromosome_6_expanded <- log_calculations(chromosome_6_expanded)
chromosome_7_expanded <- log_calculations(chromosome_7_expanded)
chromosome_8_expanded <- log_calculations(chromosome_8_expanded)
chromosome_9_expanded <- log_calculations(chromosome_9_expanded)
chromosome_10_expanded <- log_calculations(chromosome_10_expanded)
chromosome_11_expanded <- log_calculations(chromosome_11_expanded)
chromosome_12_expanded <- log_calculations(chromosome_12_expanded)
chromosome_13_expanded <- log_calculations(chromosome_13_expanded)
chromosome_14_expanded <- log_calculations(chromosome_14_expanded)
chromosome_15_expanded <- log_calculations(chromosome_15_expanded)
chromosome_16_expanded <- log_calculations(chromosome_16_expanded)
chromosome_17_expanded <- log_calculations(chromosome_17_expanded)


#List whole genome

genome_sample <- list(chromosome_1_expanded,
                      chromosome_2_expanded,
                      chromosome_3_expanded,
                      chromosome_4_expanded,
                      chromosome_5_expanded,
                      chromosome_6_expanded,
                      chromosome_7_expanded,
                      chromosome_8_expanded,
                      chromosome_9_expanded,
                      chromosome_10_expanded,
                      chromosome_11_expanded,
                      chromosome_12_expanded,
                      chromosome_13_expanded,
                      chromosome_14_expanded,
                      chromosome_15_expanded,
                      chromosome_16_expanded,
                      chromosome_17_expanded)

#Whole genome as data frame
whole_genome <- rbindlist(genome_sample) %>% as.data.frame()


#Sample whole genome
genome_sample <- lapply(genome_sample, 
                        function(x) {
                          s <- c(seq(from = 1, to = nrow(x), by = 300))
                          data <- x[c(s),]
                          return(data)
                        })



#Applying RA1
#250000
ch1 <- RA1(chromosome_1_expanded, "FPKM_250000", 250000 )
ch2 <- RA1(chromosome_2_expanded, "FPKM_250000", 250000 )
ch3 <- RA1(chromosome_3_expanded, "FPKM_250000", 250000 )
ch4 <- RA1(chromosome_4_expanded, "FPKM_250000", 250000 )
ch5 <- RA1(chromosome_5_expanded, "FPKM_250000", 250000 )
ch6 <- RA1(chromosome_6_expanded, "FPKM_250000", 250000 )
ch7 <- RA1(chromosome_7_expanded, "FPKM_250000", 250000 )
ch8 <- RA1(chromosome_8_expanded, "FPKM_250000", 250000 )
ch9 <- RA1(chromosome_9_expanded, "FPKM_250000", 250000 )
ch10 <- RA1(chromosome_10_expanded, "FPKM_250000", 250000 )
ch11 <- RA1(chromosome_11_expanded, "FPKM_250000", 250000 )
ch12 <- RA1(chromosome_12_expanded, "FPKM_250000", 250000 )
ch13 <- RA1(chromosome_13_expanded, "FPKM_250000", 250000 )
ch14 <- RA1(chromosome_14_expanded, "FPKM_250000", 250000 )
ch15 <- RA1(chromosome_15_expanded, "FPKM_250000", 250000 )
ch16 <- RA1(chromosome_16_expanded, "FPKM_250000", 250000 )
ch17 <- RA1(chromosome_17_expanded, "FPKM_250000", 250000 )


#Applying RA2
#250000
ch1_2 <- RA2(chromosome_1_expanded, "FPKM_250000_2", 250000 )
ch2_2 <- RA2(chromosome_2_expanded, "FPKM_250000_2", 250000 )
ch3_2 <- RA2(chromosome_3_expanded, "FPKM_250000_2", 250000 )
ch4_2 <- RA2(chromosome_4_expanded, "FPKM_250000_2", 250000 )
ch5_2 <- RA2(chromosome_5_expanded, "FPKM_250000_2", 250000 )
ch6_2 <- RA2(chromosome_6_expanded, "FPKM_250000_2", 250000 )
ch7_2 <- RA2(chromosome_7_expanded, "FPKM_250000_2", 250000 )
ch8_2 <- RA2(chromosome_8_expanded, "FPKM_250000_2", 250000 )
ch9_2 <- RA2(chromosome_9_expanded, "FPKM_250000_2", 250000 )
ch10_2 <- RA2(chromosome_10_expanded, "FPKM_250000_2", 250000 )
ch11_2 <- RA2(chromosome_11_expanded, "FPKM_250000_2", 250000 )
ch12_2 <- RA2(chromosome_12_expanded, "FPKM_250000_2", 250000 )
ch13_2 <- RA2(chromosome_13_expanded, "FPKM_250000_2", 250000 )
ch14_2 <- RA2(chromosome_14_expanded, "FPKM_250000_2", 250000 )
ch15_2 <- RA2(chromosome_15_expanded, "FPKM_250000_2", 250000 )
ch16_2 <- RA2(chromosome_16_expanded, "FPKM_250000_2", 250000 )
ch17_2 <- RA2(chromosome_17_expanded, "FPKM_250000_2", 250000 )




#Applying RA1
#100000
ch1 <- RA1(ch1, "FPKM_100000", 100000 )
ch2 <- RA1(ch2, "FPKM_100000", 100000 )
ch3 <- RA1(ch3, "FPKM_100000", 100000 )
ch4 <- RA1(ch4, "FPKM_100000", 100000 )
ch5 <- RA1(ch5, "FPKM_100000", 100000 )
ch6 <- RA1(ch6, "FPKM_100000", 100000 )
ch7 <- RA1(ch7, "FPKM_100000", 100000 )
ch8 <- RA1(ch8, "FPKM_100000", 100000 )
ch9 <- RA1(ch9, "FPKM_100000", 100000 )
ch10 <- RA1(ch10, "FPKM_100000", 100000 )
ch11 <- RA1(ch11, "FPKM_100000", 100000 )
ch12 <- RA1(ch12, "FPKM_100000", 100000 )
ch13 <- RA1(ch13, "FPKM_100000", 100000 )
ch14 <- RA1(ch14, "FPKM_100000", 100000 )
ch15 <- RA1(ch15, "FPKM_100000", 100000 )
ch16 <- RA1(ch16, "FPKM_100000", 100000 )
ch17 <- RA1(ch17, "FPKM_100000", 100000 )


#Applying RA2
#100000
ch1_2 <- RA2(ch1_2, "FPKM_100000_2", 100000 )
ch2_2 <- RA2(ch2_2, "FPKM_100000_2", 100000 )
ch3_2 <- RA2(ch3_2, "FPKM_100000_2", 100000 )
ch4_2 <- RA2(ch4_2, "FPKM_100000_2", 100000 )
ch5_2 <- RA2(ch5_2, "FPKM_100000_2", 100000 )
ch6_2 <- RA2(ch6_2, "FPKM_100000_2", 100000 )
ch7_2 <- RA2(ch7_2, "FPKM_100000_2", 100000 )
ch8_2 <- RA2(ch8_2, "FPKM_100000_2", 100000 )
ch9_2 <- RA2(ch9_2, "FPKM_100000_2", 100000 )
ch10_2 <- RA2(ch10_2, "FPKM_100000_2", 100000 )
ch11_2 <- RA2(ch11_2, "FPKM_100000_2", 100000 )
ch12_2 <- RA2(ch12_2, "FPKM_100000_2", 100000 )
ch13_2 <- RA2(ch13_2, "FPKM_100000_2", 100000 )
ch14_2 <- RA2(ch14_2, "FPKM_100000_2", 100000 )
ch15_2 <- RA2(ch15_2, "FPKM_100000_2", 100000 )
ch16_2 <- RA2(ch16_2, "FPKM_100000_2", 100000 )
ch17_2 <- RA2(ch17_2, "FPKM_100000_2", 100000 )



#Applying RA1
#50000
ch1 <- RA1(ch1, "FPKM_50000", 50000 )
ch2 <- RA1(ch2, "FPKM_50000", 50000 )
ch3 <- RA1(ch3, "FPKM_50000", 50000 )
ch4 <- RA1(ch4, "FPKM_50000", 50000 )
ch5 <- RA1(ch5, "FPKM_50000", 50000 )
ch6 <- RA1(ch6, "FPKM_50000", 50000 )
ch7 <- RA1(ch7, "FPKM_50000", 50000 )
ch8 <- RA1(ch8, "FPKM_50000", 50000 )
ch9 <- RA1(ch9, "FPKM_50000", 50000 )
ch10 <- RA1(ch10, "FPKM_50000", 50000 )
ch11 <- RA1(ch11, "FPKM_50000", 50000 )
ch12 <- RA1(ch12, "FPKM_50000", 50000 )
ch13 <- RA1(ch13, "FPKM_50000", 50000 )
ch14 <- RA1(ch14, "FPKM_50000", 50000 )
ch15 <- RA1(ch15, "FPKM_50000", 50000 )
ch16 <- RA1(ch16, "FPKM_50000", 50000 )
ch17 <- RA1(ch17, "FPKM_50000", 50000 )


#Applying RA2
#50000
ch1_2 <- RA2(ch1_2, "FPKM_50000_2", 50000 )
ch2_2 <- RA2(ch2_2, "FPKM_50000_2", 50000 )
ch3_2 <- RA2(ch3_2, "FPKM_50000_2", 50000 )
ch4_2 <- RA2(ch4_2, "FPKM_50000_2", 50000 )
ch5_2 <- RA2(ch5_2, "FPKM_50000_2", 50000 )
ch6_2 <- RA2(ch6_2, "FPKM_50000_2", 50000 )
ch7_2 <- RA2(ch7_2, "FPKM_50000_2", 50000 )
ch8_2 <- RA2(ch8_2, "FPKM_50000_2", 50000 )
ch9_2 <- RA2(ch9_2, "FPKM_50000_2", 50000 )
ch10_2 <- RA2(ch10_2, "FPKM_50000_2", 50000 )
ch11_2 <- RA2(ch11_2, "FPKM_50000_2", 50000 )
ch12_2 <- RA2(ch12_2, "FPKM_50000_2", 50000 )
ch13_2 <- RA2(ch13_2, "FPKM_50000_2", 50000 )
ch14_2 <- RA2(ch14_2, "FPKM_50000_2", 50000 )
ch15_2 <- RA2(ch15_2, "FPKM_50000_2", 50000 )
ch16_2 <- RA2(ch16_2, "FPKM_50000_2", 50000 )
ch17_2 <- RA2(ch17_2, "FPKM_50000_2", 50000 )



#List full genome of with MAs

genome_1_sample <- list(ch1, 
                        ch2, 
                        ch3, 
                        ch4, 
                        ch5, 
                        ch6, 
                        ch7, 
                        ch8,
                        ch9,
                        ch10,
                        ch11,
                        ch12,
                        ch13,
                        ch14,
                        ch15,
                        ch16,
                        ch17) 



genome_2_sample <- list(ch1_2, 
                        ch2_2, 
                        ch3_2, 
                        ch4_2, 
                        ch5_2, 
                        ch6_2, 
                        ch7_2, 
                        ch8_2,
                        ch9_2,
                        ch10_2,
                        ch11_2,
                        ch12_2,
                        ch13_2,
                        ch14_2,
                        ch15_2,
                        ch16_2,
                        ch17_2) 


#Sample RA1 full genome
genome_1_sample <- lapply(genome_1_sample, 
                          function(x) {
                            s <- c(seq.int(from = 1, to = nrow(x), by = 300))
                            x[c(s),]
                          })
#Sample RA2 full genome
genome_2_sample <- lapply(genome_2_sample, 
                          function(x) {
                            s <- c(seq.int(from = 1, to = nrow(x), by = 300))
                            x[c(s),]
                          })






