setwd("~/RNA/Daniel/")

load(file = "FcR_20April.RData")

#Set up

Fc_interest <- Fc_exprs3 %>% filter(Fc_Receptor != "FCRL2", Fc_Receptor != "FCRL3", Fc_Receptor != "FCRL4", Fc_Receptor != "FCRL5", Fc_Receptor != "FCRL6", Fc_Receptor != "FCGR2A.1", Fc_Receptor != "FCGBP", Fc_Receptor != "FCGRT")

Fc_interest$Diagnosis <- as.factor(Fc_interest$Diagnosis)
levels(Fc_interest2$Diagnosis) <- c("nTD", "TD")

Fc_interest2 <- Fc_interest %>% filter(Time != -32, Time != -30,  Time != 180, Time !=2, Time !=3, Time !=4, Time !="NA")
levels(Fc_interest2$Diagnosis) <- c("nTD", "TD")


#Add in t-test/wilcoxon

my_comparisons <- list(c("nTD","TD"))

Fc_interest2 %>% filter(Time == 9) %>%
  ggplot(aes(x= Diagnosis, y= Exprs, group = Diagnosis, fill = Diagnosis)) + geom_boxplot()+ facet_wrap(~Fc_Receptor) +
  labs(y = "Normalised Gene Expression", x = "Diagnosis", title = "Day 9") +theme_bw() +scale_fill_manual(values=c("seagreen3", "sandybrown")) + theme(legend.position = "none", axis.text.x = element_text(face="bold")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", vjust = 0.45, hide.ns = TRUE)

?stat_compare_means
#Baseline
Fc_interest2 %>% filter(Time == 1) %>%
  ggplot(aes(x= Diagnosis, y= Exprs, group = Diagnosis, fill = Diagnosis)) + geom_boxplot()+ facet_wrap(~Fc_Receptor) +
  labs(y = "Normalised Gene Expression", x = "Diagnosis", title = "Day 1") +theme_bw() +scale_fill_manual(values=c("seagreen3", "sandybrown")) + theme(legend.position = "none", axis.text.x = element_text(face="bold")) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format", vjust = 0.25, hide.ns = TRUE)

#Set up new loop
Fc_interest2$Time <- as.factor(Fc_interest2$Time)
levels(Fc_interest2$Time)
time_point <- unique(Fc_interest2$Time)
time_plots <- list()

for(i in levels(Fc_interest2$Time)){
  time_plots[[i]] = ggplot(Fc_interest2 %>% filter(Time == i), aes(x=Diagnosis, y= Exprs, group = Diagnosis, fill = Diagnosis)) + geom_boxplot() + facet_wrap(~Fc_Receptor) +theme_bw() +scale_fill_manual(values=c("seagreen3", "sandybrown")) + theme(legend.position = "none", axis.text.x = element_text(face="bold")) + labs(y = "Normalised Gene Expression", x = "Diagnosis", title =paste("Day",i, sep = " "))  +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", vjust = 0.45, hide.ns = TRUE)
  print(time_plots[[i]])
  Sys.sleep(1)
}

#Overall changes ####
#test
Fc_interest2 %>% filter(Diagnosis != "NA", Fc_Receptor == "FCGR1B") %>% 
  ggplot(aes(x= Time, y= Exprs, group = Time, colour = Diagnosis)) + geom_boxplot() + geom_point() + scale_x_continuous(limits = c(0, 16),breaks = seq(0,28,2))+ labs(y = "Normalised Gene Expression", x = "Days Post-Challenge", title = "FCGR1B") +scale_colour_manual(values=c("seagreen3", "sandybrown"))

#all FcRs across time
Fc_interest2 %>% filter(Diagnosis != "NA") %>% 
  ggplot(aes(x= Time, y= Exprs, group = Time, colour = Diagnosis)) + geom_boxplot(outlier.shape = NA) + geom_point(size = 1.4) + scale_x_continuous(limits = c(0, 16), breaks = seq(0,28,2)) + facet_wrap(~Fc_Receptor) +
  labs(y = "Normalised Gene Expression", x = "Days Post-Challenge") +scale_colour_manual(values=c("seagreen3", "sandybrown"))



#Separate graphs for TD participants and nTD
Fc_interest2 %>% filter(Diagnosis == "TD") %>% 
  ggplot(aes(x= Time, y= Exprs, group = Time, colour = Diagnosis)) + geom_boxplot(outlier.shape = NA) + scale_x_continuous(limits = c(-0.5, 16.5), breaks = seq(0,28,2)) + facet_wrap(~Fc_Receptor) +
  labs(y = "Normalised Gene Expression", x = "Days Post-Challenge", title = "FcR Expression in TD Participants") +scale_colour_manual(values=c("sandybrown")) + theme(legend.position = "none")

Fc_interest2 %>% filter(Diagnosis == "nTD") %>% 
  ggplot(aes(x= Time, y= Exprs, group = Time, colour = Diagnosis)) + geom_boxplot(outlier.shape = NA) + scale_x_continuous(limits = c(-0.5, 16.5), breaks = seq(0,28,2)) + facet_wrap(~Fc_Receptor) +
  labs(y = "Normalised Gene Expression", x = "Days Post-Challenge", title = "FcR Expression in nTD Participants") +scale_colour_manual(values=c("seagreen3")) + theme(legend.position = "none")


#Day of diagnosis highlight
Fc_interest2 %>% filter(Diagnosis != "NA") %>%
  mutate(TD = (Timepoint3 == "TD") )%>%
  ggplot(aes(x= Time, y= Exprs, group = Time, colour = TD)) + geom_boxplot(outlier.shape = NA) + geom_point(size =1.4) + scale_x_continuous(limits = c(-0.5, 16.5)) + facet_wrap(~Fc_Receptor) +
  labs(y = "Normalised Expression", x = "Days Post-Challenge") +scale_colour_manual(values=c("seagreen3", "sandybrown"))

##Plot sizes = 1300/600 #
save.image(file = "FcR_20April.RData")

#Plot FcGR2a against FCRLB expression ####

#Reformat to back to wide

df_t <- transpose(symbol_exprs)
gnames <- df_t[1,] %>% as.vector(mode = "character")
names(df_t) <- gnames
df_t <- df_t[-c(1),] 
Fc_exprs <-sapply(df_t, as.numeric)
pheno_exprs <- cbind(pData, Fc_exprs)

Fc_interest <- pheno_exprs %>% dplyr::select(c("Part_number","Time","FCGR3A","FCAR","FCGR1B","FCRLA","FCRLB","FCGR2A")) %>%
  filter(Time == "7")#89 row
#D7
Fc_interest %>% ggplot(aes(x=FCGR2A,y=FCRLB)) + geom_point()  +
  geom_smooth(method=lm) +
  stat_cor(method = "spearman") + theme_bw() +labs(title = "FcR expression at D7")

