#This should be run after applying the Setup and Function code


#Setting theme
theme_set(theme_classic(base_size = 14,
                        base_family = 'Arial',
                        base_line_size = 0.5,
                        base_rect_size = 0.5))


brewer.pal(6,'Dark2')
display.brewer.pal(6,'Dark2')





##FINDING TOP_N GENES IN GENOME

#Top 10 highest constitutively expressing genes by mean

top_10_mean <- top_n(data_frame_2, 10, mean_FPKM)
top_10_mean_log <- top_n(data_frame_2, 10, log_FPKM)


# Top 10 highest constitutively expressing genes by sum

top_10_sum <- top_n(data_frame_2, 10, sum_FPKM)
top_10_sum_log <- top_n(data_frame_2, 10, log_sum)


# Top 10 lowest variance genes

bottom_10_var <- top_n(data_frame_2, -10, CV)
bottom_10_var_log <- top_n(data_frame_2, -10, log_var)






#PLOTTING VARIATION ~ FPKM

#log_var ~ log_FPKM
var_FPKM_log <- ggplot(genome_sample, aes(x = log_FPKM, y = log_var)) +
  geom_point(colour = "#D4D4D4" ) +
  geom_point(data = top_10_mean_log, aes(x = log_FPKM, y = log_var), colour = "#FBB962") +
  geom_point(data = bottom_10_var_log, aes(x = log_FPKM, y = log_var), colour = "#51CDC9") +
  geom_smooth(method = 'lm', se = FALSE, colour = "black", size = 0.5) +
  xlab("log10(mean FPKM)") + ylab("log10(variance)") +
  theme(axis.title = element_text(size = 10))

ggsave(filename = "log(var)_log(mean)_candidates.png", var_FPKM_log, width = 178, height = 178*0.6, unit = "mm")

#log_sd ~ log_FPKM
sd_FPKM_log <- ggplot(genome_sample, aes(x = log_FPKM, y = log_sd)) +
  geom_point(colour = "#D4D4D4" ) +
  geom_point(data = top_10_mean_log, aes(x = log_FPKM, y = log_var), colour = "#FBB962") +
  geom_point(data = bottom_10_var_log, aes(x = log_FPKM, y = log_var), colour = "#51CDC9") +
  geom_smooth(method = 'lm', se = FALSE, colour = "black", size = 0.5) +
  xlab("log10(mean FPKM)") + ylab("log10(sd)") +
  theme(axis.title = element_text(size = 10))

ggsave(filename = "log(sd)_log(mean)_candidates.png", sd_FPKM_log, width = 178, height = 178*0.6, unit = "mm")

#CV ~ mean_FPKM
var_FPKM <- ggplot(data = genome_sample, aes(x = mean_FPKM, y = CV )) +
  geom_point(colour = "#D4D4D4") +
  geom_point(data = top_10_mean, colour = "#FBB962") +
  geom_point(data = bottom_10_var, colour = "#51CDC9") +
  xlab("mean FPKM") + ylab("CV") +
  theme(axis.title = element_text(size = 10)) +
  coord_cartesian(expand = FALSE) +
  xlim(c(0, 8000))

ggsave(filename = "CV_log(mean)_candidates.png", var_FPKM, width = 178, height = 178*0.6, unit = "mm")




#PEARSONS CORRELATION COEFFICIENT FOR log_var ~ log_mean
cor.test(whole_genome$log_FPKM, whole_genome$log_var, method = "pearson")


#MA PLOTS 250000

#Genes previously targeted with high transformation efficiency
high <- c("Cre03.g161400.v5.5", "Cre13.g586300.v5.5", "Cre05.g241450.v5.5")

#Genes previously targeted with low transformation efficiency
low <- c("Cre02.g082550.v5.5", "Cre12.g498550.v5.5")

#Plot
s <- seq(from = -7.5, to = 4.5, by = 0.5)

ma_plot <- function(data_1, data_2, data){
  
  a <- data_1
  b <- data_2
  c <- data
  
  d <- ggplot(c, mapping = aes(x = bp_position), colour = "#D4D4D4") +
    geom_point(mapping = aes(y = log_FPKM), colour = "#D4D4D4") + 
    geom_point(data = c[is.na(c$gene_short_name),],  mapping = aes(y = log_FPKM), colour = "darkgray") +
    geom_point(data = c[c$gene_id %in%  high,], mapping = aes(y = log_FPKM), colour = "#DD0000") +
    geom_point(data = c[c$gene_id %in%  low,], mapping = aes(y = log_FPKM), colour = "#296d98" ) +
    geom_line(data = b, mapping = aes(y = FPKM_250000_2), na.rm = TRUE, colour = "#D95F02" ) +
    geom_line(data = a, mapping = aes(y = FPKM_250000), na.rm = TRUE, colour = "#1B9E77" ) +
    facet_wrap(~chromosome) +
    scale_y_continuous(breaks = s) +
    labs(title = "FPKM_250000")
  
  
  return(d)
}

x <- Map(f = ma_plot, data_1 = genome_1_sample, data_2 = genome_2_sample, data = genome_sample)

for(i in 1:17){
  ggsave(filename = paste("FPKM_250000",i,".png", sep = " "), x[[i]], width = 178, height = 178*0.6, unit = "mm")
}



#MA PLOTS 100000
s <- seq(from = -7.5, to = 4.5, by = 0.5)

ma_plot <- function(data_1, data_2, data, llim, ulim){
  
  a <- data_1
  b <- data_2
  c <- data
  
  d <- ggplot(c, mapping = aes(x = bp_position), colour = "#D4D4D4") +
    geom_point(mapping = aes(y = log_FPKM), colour = "#D4D4D4") +
    geom_point(data = c[is.na(c$gene_short_name),],  mapping = aes(y = log_FPKM), colour = "darkgray") +
    geom_point(data = c[c$gene_id %in% high,], mapping = aes(y = log_FPKM), colour = "#DD0000") +
    geom_point(data = c[c$gene_id %in% low,], mapping = aes(y = log_FPKM), colour = "#296d98" ) +
    geom_line(data = b, mapping = aes(y = FPKM_100000_2), na.rm = TRUE, colour = "#E7298A" ) +
    geom_line(data = a, mapping = aes(y = FPKM_100000), na.rm = TRUE, colour = "#7570B3" ) +
    facet_wrap(~chromosome) +
    coord_cartesian(xlim = c(llim, ulim)) +
    geom_hline(yintercept = 0.5, colour = "#E7298A", linetype = "dashed") +
    scale_y_continuous(breaks = s) +
    labs(title = "FPKM_100000") + ylab("log10(FPKM)") + xlab("bp position along chromosome")
  
  
  return(d)
}


f <- Map(f = ma_plot, data_1 = genome_1_sample, data_2 = genome_2_sample, data = genome_sample)

for(i in 1:17){
  ggsave(filename = paste("FPKM_100000",i,".png", sep = " "), f[[i]], width = 178, height = 178*0.6, unit = "mm")
}




#MA PLOTS 50000
s <- seq(from = -7.5, to = 4.5, by = 0.5)

ma_plot <- function(data_1, data_2, data, llim, ulim){
  
  a <- data_1
  b <- data_2
  c <- data
  d <- a[(a$log_FPKM > 2), ] 
  d <- d[!duplicated(d$gene_short_name),]
  
  e <- ggplot(c, mapping = aes(x = bp_position), colour = "#D4D4D4") +
    geom_point(mapping = aes(y = log_FPKM), colour = "#D4D4D4") +
    geom_point(data = c[is.na(c$gene_short_name),],  mapping = aes(y = log_FPKM), colour = "darkgray") +
    geom_point(data = c[c$gene_id %in% high,], mapping = aes(y = log_FPKM), colour = "#DD0000") +
    geom_point(data = c[c$gene_id %in% low,], mapping = aes(y = log_FPKM), colour = "#296d98" ) +
    geom_line(data = b, mapping = aes(y = FPKM_50000_2), na.rm = TRUE, colour = "#E6AB02" ) +
    geom_line(data = a, mapping = aes(y = FPKM_50000), na.rm = TRUE, colour = "#66A61E" ) +
    facet_wrap(~chromosome) +
    coord_cartesian(xlim = c(llim, ulim)) +
    geom_hline(yintercept = 0.5, colour = "#E7298A", linetype = "dashed") +
    geom_text(data = d, aes(x = bp_position, y = log_FPKM, label = gene_short_name), check_overlap = TRUE, size = 2.5, nudge_y = -0.1) +
    scale_y_continuous(breaks = s) +
    labs(title = "Moving Average Plot with 50000 Rolling Window") + xlab("bp along chromosome") + ylab("log10(FPKM)")
  
  return(e)
}


j <- Map(f = ma_plot, data_1 = genome_1_sample, data_2 = genome_2_sample, data = genome_sample)

for(i in 1:17){
  ggsave(filename = paste("FPKM_50000",i,".png", sep = " "), j[[i]], width = 178, height = 178*0.6, unit = "mm")
}





#CALCULATING AREA UNDER CURVE

#Finding Intercepts Function
intercepts <- function(i, FPKM, data) {
  df <- data
  l <- ifelse((df[i, FPKM] <= 0 && df[i + 1, FPKM] >= 0) || (df[i, FPKM] >= 0 && df[i + 1, FPKM] <= 0), df[i, 'bp_position'], FALSE)
  return(l)
}




#FUNCTION TO FIND Y = 0 INTERCEPTS

intercepts <- function(i, FPKM, data) {
  df <- data
  l <- ifelse((df[i, FPKM] <= 0 && df[i + 1, FPKM] >= 0) || (df[i, FPKM] >= 0 && df[i + 1, FPKM] <= 0), df[i, 'bp_position'], FALSE)
  return(l)
}


#Applying intercepts function for 250000
intersection_points <- lapply(seq_along(genome_1_sample), function(x){
  p <- sapply(1:(nrow(genome_1_sample[[x]]) - 1), intercepts, 20, genome_1_sample[[x]]) %>% unlist()
  v <- p[p > 0]
  v <- na.omit(v)
  return(as.vector(v))
} )


#Subseting data for FPKM_100000

subbed_100 <- list()

subbed_100[[1]] <- limit(genome_1_sample[[1]], 2292765, 5771765)
subbed_100[[2]] <- limit(genome_1_sample[[2]], 4039905, 7514905)
subbed_100[[3]] <- limit(genome_1_sample[[3]], 4233119, 6102119)
subbed_100[[4]] <- limit(genome_1_sample[[6]], 1188148, 3382648)
subbed_100[[5]] <- limit(genome_1_sample[[9]], 2536233, 4430733)
subbed_100[[6]] <- limit(genome_1_sample[[10]],  1218313, 4199313)
subbed_100[[7]] <- limit(genome_1_sample[[12]],  2814765, 4840265)
subbed_100[[8]] <- limit(genome_1_sample[[12]],  5031765, 6663765)
subbed_100[[9]] <- limit(genome_1_sample[[16]],  346445, 3320445)
subbed_100[[10]] <- limit(genome_1_sample[[17]], 2637754, 4125254)


#Finding intersection points for 100000
intersection_points <- lapply(seq_along(genome_1_sample), function(x){
  p <- sapply(1:(nrow(genome_1_sample[[x]]) - 1), intercepts, 20, genome_1_sample[[x]]) %>% unlist()
  v <- p[p > 0]
  v <- na.omit(v)
  return(as.vector(v))
} )


#Join intercepts to data function
join <- function (data_1, data_2) {
  c <- data.frame(bp_position = data_1,
                  intercept = data_1)
  d <- left_join(c, as.data.frame(data_2), by = "bp_position")
  d <- subset(d, select = c('bp_position', 'intercept'))
  d <- left_join(data_2, d, by = "bp_position")
  d$area <- NA
  return(d)
}

#Applying join function - Match intersection bp to FPKM_250000
t <- Map(f = join, data_1 = intersection_points, data_2 = genome_1_sample)

t <- Map(f = join, data_1 = intersection_points, data_2 = subbed_100)



#Area under curve function

area <- function(i, data){
  df <- data
  r <- data[, 'intercept'] %>% na.omit() %>% as.vector()
  
  a <- MESS::auc(x = df[,'bp_position'], y = df[,'FPKM_250000'], from = r[i], to = r[i + 1])
  
  data[data$bp_position >= r[i] & data$bp_position <= r[i + 1], 'area'] <- a
  
  return(data)
}

#Applying area function

area_func <- function(x){
  r <- t[[x]][, 'intercept'] %>% na.omit() %>% as.vector()
  
  for(i in 1:(length(r)-1)){
    t[[x]] <- area(i, t[[x]])
  }

  return(t[[x]])
}


#Appply area function
g <- lapply(seq_along(t), area_func)


#Finding the largest areas

largest_area <- lapply(seq_along(g), function(x){
  a <- g[[x]][, c('area', 'bp_position', 'chromosome', 'intercept', 'gene_id', 'gene_short_name', 'FPKM_250000')]
  a <- a %>% drop_na("area")
  a <- a[a$area > 0, ]
  return(a)
})


la <- largest_area %>% rbindlist()

la <- la[!duplicated(la[,c('area')]),]


#Top 10 area 
top_10_area <- la %>% top_n(10, area) 

#250000
top_10_area$end_intercept <- c(7514905, 5771765, 3320445, 4199313, 4840265, 3382648, 4430733, 6663765, 6102119, 4125254)


#100000
top_10_area$end_intercept <- c(3152265 , 5568405, 7495405, 5569119, 3291648, 3845233, 3787765, 4717265, 3269445, 3886254)

#Finding the largest areas

largest_area <- lapply(seq_along(g), function(x){
  a <- g[[x]][, c('area', 'bp_position', 'chromosome', 'intercept', 'gene_id', 'gene_short_name', 'FPKM_250000')]
  a <- a %>% drop_na("area")
  a <- a[a$area > 0, ]
  return(a)
})


la <- largest_area %>% rbindlist()

la <- la[!duplicated(la[,c('area')]),]


#Top 10 area 
top_10_area <- la %>% top_n(10, area) 

#Export to Excel
write_xlsx(top_10_area,"/Users/traceyly/Desktop/R Project/Time Data\\top_10_area 100000.xlsx")



#ZOOMED IN PLOTS

#100000
f = list()

f[[1]] <- ma_plot(genome_1_sample[[1]], genome_2_sample[[1]], genome_sample[[1]], 2292765, 5771765)
f[[2]] <- ma_plot(genome_1_sample[[2]], genome_2_sample[[2]], genome_sample[[2]], 4039905, 7514905)
f[[3]] <- ma_plot(genome_1_sample[[3]], genome_2_sample[[3]], genome_sample[[3]], 4233119, 6102119)
f[[4]] <- ma_plot(genome_1_sample[[6]], genome_2_sample[[6]], genome_sample[[6]], 1188148, 3382648)
f[[5]] <- ma_plot(genome_1_sample[[9]], genome_2_sample[[9]], genome_sample[[9]], 2536233, 4430733)
f[[6]] <- ma_plot(genome_1_sample[[10]], genome_2_sample[[10]],genome_sample[[10]], 1218313, 4199313)
f[[7]] <- ma_plot(genome_1_sample[[12]], genome_2_sample[[12]],genome_sample[[12]], 2814765, 4840265)
f[[8]] <- ma_plot(genome_1_sample[[12]], genome_2_sample[[12]],genome_sample[[12]], 5031765, 6663765)
f[[9]] <- ma_plot(genome_1_sample[[16]], genome_2_sample[[16]],genome_sample[[16]], 346445, 3320445)
f[[10]] <- ma_plot(genome_1_sample[[17]], genome_2_sample[[17]], genome_sample[[17]], 2637754, 4125254)

for(i in 1:10){
  ggsave(filename = paste("FPKM_100000_zoomed",i,".png", sep = " "), f[[i]], width = 178, height = 178*0.6, unit = "mm")
}

#50000

h = list()

h[[1]] <- ma_plot(genome_1_sample[[1]], genome_2_sample[[1]], genome_sample[[1]], 2367265, 3152265) 
h[[2]] <- ma_plot(genome_1_sample[[2]], genome_2_sample[[2]], genome_sample[[2]], 4605905, 5568405)
h[[3]] <- ma_plot(genome_1_sample[[2]], genome_2_sample[[2]], genome_sample[[2]], 6632905, 7495405)
h[[4]] <- ma_plot(genome_1_sample[[3]], genome_2_sample[[3]], genome_sample[[3]], 4876119, 5569119)
h[[5]] <- ma_plot(genome_1_sample[[6]], genome_2_sample[[6]], genome_sample[[6]], 1660148, 3291648) 
h[[6]] <- ma_plot(genome_1_sample[[9]], genome_2_sample[[9]], genome_sample[[9]], 2755233, 3845233) 
h[[7]] <- ma_plot(genome_1_sample[[12]], genome_2_sample[[12]], genome_sample[[12]], 2928765, 3787765) 
h[[8]] <- ma_plot(genome_1_sample[[12]], genome_2_sample[[12]], genome_sample[[12]], 3800265, 4717265) 
h[[9]] <- ma_plot(genome_1_sample[[16]], genome_2_sample[[16]], genome_sample[[16]], 2713445, 3269445) 
h[[10]] <- ma_plot(genome_1_sample[[17]], genome_2_sample[[17]], genome_sample[[17]], 2932254, 3886254) 


for(i in 1:10){
  ggsave(filename = paste("FPKM_50000_zoomed",i,".png", sep = " "), h[[i]], width = 178, height = 178*0.6, unit = "mm")
}

#PRELIMINARY CANDIDATES

#Preliminary candidates
prelim_candidates <- c("Cre01.g013450",
                       "Cre01.g013700",
                       "RBCS2",
                       "RBCS1",
                       "Cre09.g388726",
                       "ASD1",
                       "RPL24",
                       "RPS7",
                       "RPS17",
                       "Cre12.g498600",
                       "Cre12.g499352",
                       "PGH1",
                       "RPS11",
                       "APR1",
                       "RPP0",
                       "PRPS6",
                       "Cre12.g522600",
                       "Cre16.g663900",
                       "ASA5"
)

#Table of preliminary candidates
prelim_candidates <- data_frame[data_frame$gene_short_name %in% prelim_candidates,]


#TIME DISTRIBUTION PLOTS

#Candidates time distribution

f <- c("dark.11h", "dark.9h", "dark.7h", "dark.5h", "dark.3h", "dark.1h", "dark.0.5h", "light.0h", "light.0.5h", "light.1h", "light.3h", "light.5h", "light.7h", "light.9h", "light.11h", "dark.13h" )

candidates <- c("Cre02.g120100.v5.5", "Cre02.g120150.v5.5", "Cre09.g388726.v5.5", "Cre12.g499352.v5.5", "Cre12.g517150.v5.5" )

candidates <- data_frame[data_frame$gene_id %in% candidates,]

candidates <- candidates %>% gather(f, key = "time", value = "FPKM")

candidates$time <- factor(candidates$time, levels = f)

candidate_plot <- ggplot(candidates, aes(x = time, y = FPKM, colour = gene_short_name, group = gene_short_name)) +
  geom_line() +
  geom_point() +
  xlab("Time") + ylab("FPKM") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 10))

#Export table to Excel
write.table(x = candidates, file = "Candidates.xls")


#Setting up data frame

candidates <- c("Cre02.g120100.v5.5", "Cre02.g120150.v5.5", "Cre09.g388726.v5.5", "Cre12.g499352.v5.5", "Cre12.g517150.v5.5" )

targeted <- c("Cre13.g586300.v5.5", "Cre05.g241450.v5.5", "Cre03.g161400.v5.5")

f <- c("dark.11h", "dark.9h", "dark.7h", "dark.5h", "dark.3h", "dark.1h", "dark.0.5h", "light.0h", "light.0.5h", "light.1h", "light.3h", "light.5h", "light.7h", "light.9h", "light.11h", "dark.13h" )

k <- data_frame[data_frame$gene_id %in% c(candidates, targeted), ]

k$time_var <- rowVars(k[, f])

k$mean_FPKM<- rowMeans(k[, f])

k$log_var <- k$time_var %>% log10()

k$log_FPKM<- k$mean_FPKM %>% log10()

k$max <- sapply(1:nrow(k), function(x){
  max(k[x, f])
})

k <- k %>% gather(f, key = "time", value = "FPKM")

k$time <- factor(k$time, levels = f)

k$scaled_FPKM <- (k$FPKM / k$range) * 100



#Targeted time distribution

p <- ggplot(k, aes(x = time, y = scaled_FPKM, colour = group, group = gene_short_name)) +
  geom_line() +
  geom_point() +
  xlab("Time") + ylab("FPKM") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 10))

ggsave(filename = paste("Targeted Time Distributions 7.png"), p, width = 178, height = 178*0.6, unit = "mm")




#Scaled Diel Cyle plots
scaled_targeted <- ggplot(k[k$gene_id %in% targeted,], aes(x = time, y = scaled_FPKM, colour = gene_short_name, group = gene_short_name)) +
  geom_line() +
  geom_point() +
  xlab("Time") + ylab("FPKM") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 10))


ggsave(filename = paste("Targeted Scaled.png"), scaled_targeted, width = 178, height = 178*0.6, unit = "mm")




#TARGTEED VS. CANDIDATES CV ~ mean FPKM

candidates <- c("Cre02.g120100.v5.5", "Cre02.g120150.v5.5", "Cre09.g388726.v5.5", "Cre12.g499352.v5.5", "Cre12.g517150.v5.5" )

candidates <- data_frame_2[data_frame_2$gene_id %in% candidates,]

targeted <- c("Cre13.g586300.v5.5", "Cre05.g241450.v5.5", "Cre03.g161400.v5.5")

targeted <- data_frame_2[data_frame_2$gene_id %in% targeted,]

brewer.pal(6,'Accent')
display.brewer.pal(6,'Accent')


target_candidate <- ggplot(data = genome_sample, aes(x = mean_FPKM, y = CV )) +
  geom_point(colour = "#D4D4D4") +
  geom_point(data = candidates, aes(x = mean_FPKM, y = CV), colour = "#386CB0") +
  geom_text(data = candidates, aes(x = mean_FPKM, y = CV, label = gene_short_name), check_overlap = TRUE, size = 2.5, nudge_y = 0.1, nudge_x = 400) +
  geom_point(data = targeted, aes(x = mean_FPKM, y = CV), colour = "#F0027F") +
  geom_text(data = targeted, aes(x = mean_FPKM, y = CV, label = gene_short_name), check_overlap = TRUE, size = 2.5, nudge_y = -0.1, nudge_x = 200) +
  xlab("mean FPKM") + ylab("CV") +
  theme(axis.title = element_text(size = 10)) +
  coord_cartesian(expand = FALSE) +
  xlim(c(0, 8000))


ggsave(filename = "targeted_candidates.png", target_candidate, width = 178, height = 178*0.6, unit = "mm")





