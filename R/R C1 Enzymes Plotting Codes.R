#C1 Enzymes Functions should be run before

#Setting theme
theme_set(theme_classic(base_size = 14,
                        base_family = 'Arial',
                        base_line_size = 0.5,
                        base_rect_size = 0.5))


brewer.pal(6,'Dark2')
display.brewer.pal(6,'Dark2')

#Gathering Data
b <- a

f <- c("dark.11h", "dark.9h", "dark.7h", "dark.5h", "dark.3h", "dark.1h", "dark.0.5h", "light.0h", "light.0.5h", "light.1h", "light.3h", "light.5h", "light.7h", "light.9h", "light.11h", "dark.13h" )

b <- b %>% gather(f, key = "time", value = "FPKM")



#Time distribution plots for C1 Enzymes

split_b <- split(b, b$Pathway)

my_plot <- function(data) {
  
  ggplot(data, aes(x = factor(time, levels = time_frame), y = FPKM, colour = Name, group = Name)) +
    geom_line() +
    geom_point() +
    facet_wrap(~Pathway, scales = "free") +
    labs(title = "Expression over diel cycle", x = "Time") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

p <- lapply(split_b, my_plot)

p <- ggarrange(plotlist = p, nrow = 1, ncol = 1)



#Heatmap for diel cycle

f <- c("dark.11h", "dark.9h", "dark.7h", "dark.5h", "dark.3h", "dark.1h", "dark.0.5h", "light.0h", "light.0.5h", "light.1h", "light.3h", "light.5h", "light.7h", "light.9h", "light.11h", "dark.13h" )

e_time <-  enzymes[, c("Name", "Pathway", f)] %>% as.data.frame()


z_score <- function(time, data){
  data_1 <- data
  
  
  #Calculate mean
  data_1[,"mean"] <-  sapply(1:nrow(data_1), function(x){
    m <- mean(as.numeric(data_1[x, time]))
    return(m)
  })
  
  #Calculate sd
  data_1[,"sd"] <-  sapply(1:nrow(data_1), function(x){
    s <- sd(as.numeric(data_1[x, time]))
    return(s)
  })
  
  #Calculate Z-scores
  data_1[, time] <- (data[,time] - data_1[, "mean"])/data_1[, "sd"]
  
  return(data_1)
}

e_time <- z_score(3:18, e_time)

ord <- hclust( dist(e_time[, 3:18], method = "euclidean"), method = "ward.D" )$order %>% as.vector()

e_time[, "Name"] <- factor(e_time[, "Name"], levels = e_time[ord, "Name"])

e_time_2 <- e_time %>% gather(f, key = "Time" , value = "FPKM", convert = TRUE )

e_time_2[, "Time"] <- factor(e_time_2[,"Time"], levels = f)



diel_cycle_heatmap <- ggplot(e_time_2, aes(x = Time , y = Name, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient2(low = "#D01C8B" , high = "#4DAC26", mid = "#F7F7F7", midpoint = 0, space = "Lab", name="Z-Score", na.value = "grey") +
  facet_grid(Pathway~., scales = "free", space = "free", labeller = label_wrap_gen(width=10)) +
  xlab("Time") + ylab("Genes") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill="#D4D4D4", colour = "white"),
        panel.spacing=unit(0,"cm"),
        strip.text.y = element_text(size = 6, angle = 0),
        axis.text.x = element_text(angle = 45, vjust = 1, size = 7, hjust = 1),
        axis.text.y = element_text(angle = 0, vjust = 1, size = 5, hjust = 1),
        axis.title = element_text(size = 7))

ggsave(filename = "Heatplot C1 Enzymes Diel Cycle.png", diel_cycle_heatmap, width = 178, height = 178*0.6, unit = "mm")


for(i in 1:6){
  ggsave(filename = paste("C1 Enzymes Time Distributions",i,".png", sep = ""), p[[i]], width = 178, height = 178*0.6, unit = "mm")
}


#B12 +/- boxplot 

f <- c("dark.11h", "dark.9h", "dark.7h", "dark.5h", "dark.3h", "dark.1h", "dark.0.5h", "light.0h", "light.0.5h", "light.1h", "light.3h", "light.5h", "light.7h", "light.9h", "light.11h", "dark.13h" )


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


c <- a

g <- c(ancestral_m, ancestral_p, met_m, met_p)

c <- c %>% gather(g, key = "Trial", value = "Expression")

c <- c %>% 
  separate(Trial, into = c("Condition", "Trial"), sep = -2, convert = TRUE)

c[, "Pathway"] <- factor(c[, "Pathway"], levels = c("Folate synthesis", 
                                                    "Folate cycle" , 
                                                    "Methionine cycle", 
                                                    "Poly- & Deglutamylation", 
                                                    "DNA methyltransferases",  
                                                    "RNA methyltransferases" ))




c[, "Condition"] <- factor(c[, "Condition"], levels = c("A_pB12", "M_pB12", "A_mB12", "M_mB12"))

c[,"log_Expression"] <- c[,"Expression"] %>% log10()

c[c$log_Expression == -Inf, "log_Expression"] <- NA

Boxplot <- ggplot(c, aes(x = Condition, y = log_Expression, fill = Pathway)) +
  stat_boxplot(geom ='errorbar', width = 0.5) + 
  geom_boxplot(outlier.size = 0.5) +
  geom_point(alpha = 1/2, size = 0.75) +
  xlab("Condition") + ylab("log10(FPKM)") +
  facet_grid(.~Pathway, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 8),
        panel.spacing=unit(0,"cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5)) +
  scale_fill_brewer(palette="Set2")

ggsave(filename = paste("Boxplots 2.png"), Boxplot, width = 178, height = 178*0.6, unit = "mm")


#Shapiro-Wilk normality test (sample size is < 5000)

split_c <- split(c, c$Pathway)
split_c <- split_c[c("Folate synthesis", 
                     "Folate cycle" , 
                     "Methionine cycle", 
                     "Poly- & Deglutamylation", 
                     "DNA methyltransferases",  
                     "RNA methyltransferases" )]

shapiro.test(split_c[[1]][,"log_Expression"])
shapiro.test(split_c[[2]][,"log_Expression"])
shapiro.test(split_c[[3]][,"log_Expression"])
shapiro.test(split_c[[4]][,"log_Expression"])
shapiro.test(split_c[[5]][,"log_Expression"])
shapiro.test(split_c[[6]][,"log_Expression"])

#Wilcoxon U Test

y <- c[c$Condition == c("A_mB12", "M_mB12"), ]

y[, "Condition"] <- factor(y[,"Condition"], levels = c("A_mB12", "M_mB12"))

y <- split(y, y$Pathway)
y <- y[c("Folate synthesis", 
         "Folate cycle" , 
         "Methionine cycle", 
         "Poly- & Deglutamylation", 
         "DNA methyltransferases",  
         "RNA methyltransferases" )]

wilcox.test(y[[1]][, "Expression"] ~ y[[1]][,"Condition"])
wilcox.test(y[[2]][, "Expression"] ~ y[[2]][,"Condition"])
wilcox.test(y[[3]][, "Expression"] ~ y[[3]][,"Condition"])
wilcox.test(y[[4]][, "Expression"] ~ y[[4]][,"Condition"])
wilcox.test(y[[5]][, "Expression"] ~ y[[5]][,"Condition"])
wilcox.test(y[[6]][, "Expression"] ~ y[[6]][,"Condition"])



#C1 bar graphs with error bars 

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


d <- enzymes


d_Am <- d[, c("Name", "Pathway", ancestral_m)]

d_Am$mean <- rowMeans(d_Am[,3:6])
d_Am$sd <- rowVars(d_Am[,3:6]) %>% sqrt()
d_Am$condition <- "A_mB12"
d_Am <- d_Am[, c(1, 2, 7:9 )]

d_Mm <- d[, c("Name", "Pathway", met_m)]

d_Mm$mean <- rowMeans(d_Mm[,3:6])
d_Mm$sd <- rowVars(d_Mm[,3:6]) %>% sqrt()
d_Mm$condition <- "M_mB12"
d_Mm <- d_Mm[, c(1, 2, 7:9 )]


bars <- rbind(d_Am, d_Mm)

split_d <- split(bars, bars$Pathway)

split_d <- split_d[c("Folate synthesis", 
                     "Folate cycle" , 
                     "Methionine cycle", 
                     "Poly- & Deglutamylation", 
                     "DNA methyltransferases",  
                     "RNA methyltransferases" )]

my_plot <- function(data){
  plot <- ggplot(data, aes(x=as.factor(Name), y= mean, fill= condition)) +
    geom_bar(position=position_dodge(), stat="identity", colour='black') +
    geom_errorbar(aes(ymin= mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9)) +
    facet_wrap(.~Pathway) +
    xlab("Gene") + ylab("Mean FPKM") +
    theme(strip.placement = "outside",
          axis.text.x = element_text(angle = 90, hjust=1, size = 10),
          axis.title = element_text(size = 10)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  return(plot)
}


bar <- lapply(split_d, my_plot)

for(i in 1:6){
  ggsave(filename = paste("B12 Bars",i,".png", sep = ""), bar[[i]], width = 178, height = 178*0.6, unit = "mm")
}



#Variance ~ genes plot

enzymes_var <- a


enzymes_var$Pathway <- factor(enzymes_var$Pathway, levels = c("Folate synthesis", 
                                                              "Folate cycle" , 
                                                              "Methionine cycle", 
                                                              "Poly- & Deglutamylation", 
                                                              "DNA methyltransferases",  
                                                              "RNA methyltransferases" ), ordered = TRUE) 


enzymes_var$CV <- enzymes_var$sd / enzymes_var$mean_FPKM


Variance <- ggplot(enzymes_var, aes(x = reorder(Name, CV), y = CV)) +
  geom_bar(aes(fill = Pathway), stat = "identity", position = "dodge", colour = "black") +
  xlab("C1 Enzymes") + ylab("CV") +
  facet_grid(.~Pathway, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1, size = 5),
        axis.text.y = element_text(size = 5),
        legend.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 6),
        strip.background = element_blank(),
        panel.spacing=unit(0,"cm"),
        strip.text.x = element_blank()) +
  # geom_hline(yintercept = 1.0) +
  scale_y_continuous(breaks = seq(from = -10, to = 10, by = 0.5)) +
  scale_fill_brewer(palette="Set2") +
  coord_cartesian(expand = FALSE )


ggsave(filename = "B12 Variance 2.png", Variance ,width = 178, height = 178*0.6, unit = "mm")







