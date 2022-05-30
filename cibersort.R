#Cibersort ####

#Prep mixture file ####

#Sort dups
make_clean_names(rownames(exprs2))

#Separate by time pints first

exprs_0 = exprs2[,grep("_0_", colnames(exprs2))]
exprs_1 = exprs2[,grep("_0.5_|_1_", colnames(exprs2))]
#exprs_1 = exprs2[,grep("_1_", colnames(exprs2))]
exprs_7 = exprs2[,grep("_7_", colnames(exprs2))]

#Capitalise genes
row.names(exprs_0) <- str_to_upper(row.names(exprs_0))
row.names(exprs_1) <- str_to_upper(row.names(exprs_1))
row.names(exprs_7) <- str_to_upper(row.names(exprs_7))


write.table(exprs_0, file = "T1T2_D0_cibersort.txt", sep = "\t", quote = FALSE) #in notepad add tab to header and write Gene as colname 
write.table(exprs_1, file = "T1T2_D1_cibersort.txt", sep = "\t", quote = FALSE)
write.table(exprs_7, file = "T1T2_D7_cibersort.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

check <- fread("T1T2_D0_cibersort.txt")
?write.table

#May need to remove NA genes or impute using KNN

#Results Baseline
# extra subsets ciber_0 <- fread("cibersort/CIBERSORTx_Day0/CIBERSORTxGEP_Job12_GEPs.txt")
ciber_0 <- fread("cibersort/CIBERSORTx_Day0/CIBERSORTxGEP_Job13_GEPs.txt")
gene_list <- c("FCGR3A","FCAR","FCGR1B","FCRLA","FCRLB","FCGR2A")
ciber_0fc = subset(ciber_0, GeneSymbol %in% gene_list)

#Bar Chart
ciber_0fc2 <- melt(ciber_0fc, id.vars = c("GeneSymbol"), measure.vars = c(2:11))
colnames(ciber_0fc2) <- c("FcR","Cell.Type", "Expression")

bp_D0<- ggplot(ciber_0fc2, aes(x=FcR, y=Expression, fill=Cell.Type))+
  geom_bar(width = 1, stat = "identity") + labs(title = "Baseline FcR expression", x="\nFcR") + ylim(0,120000) + theme(legend.position = "none", axis.text.x = element_text(size= 11,face="bold"))

bp_D0


#Results 12h and 24hs ####
#Shows the relative expression of FcRs on those cell types
#Daniel's data
ciber_1 <- fread("cibersort/CIBERSORTxGEP_Job9_GEPs.txt")
  gene_list <- c("FCGR3A","FCAR","FCGR1B","FCRLA","FCRLB","FCGR2A")
ciber_1fc = subset(ciber_1, GeneSymbol %in% gene_list)

#My data
ciber_1 <- fread("cibersort/CIBERSORTx_Day1/CIBERSORTxGEP_Job11_GEPs.txt")
gene_list <- c("FCGR3A","FCAR","FCGR1B","FCRLA","FCRLB","FCGR2A")
ciber_1fc = subset(ciber_1, GeneSymbol %in% gene_list)

#Bar Chart
ciber_1fc2 <- melt(ciber_1fc, id.vars = c("GeneSymbol"), measure.vars = c(2:11))
colnames(ciber_1fc2) <- c("FcR","Cell.Type", "Expression")

bp_1 <- ggplot(ciber_1fc2, aes(x=FcR, y=Expression, fill=Cell.Type))+
  geom_bar(width = 1, stat = "identity") + labs(title = "Day 1 FcR expression\n", x="\nFcR", y="") + ylim(0,120000) + theme(legend.position = "none", axis.text.x = element_text(size= 11,face="bold"))

bp_1

#Results Day 7 ####
#Could do a day7/8 combined
ciber <- fread("cibersort/CIBERSORTx_Day7/CIBERSORTxGEP_Job10_GEPs.txt")
gene_list <- c("FCGR3A","FCAR","FCGR1B","FCRLA","FCRLB","FCGR2A")
ciber_fc = subset(ciber, GeneSymbol %in% gene_list)

ciber_fc2 <- melt(ciber_fc, id.vars = c("GeneSymbol"), measure.vars = c(2:11))
colnames(ciber_fc2) <- c("FcR","Cell.Type", "Expression")

bp_7 <- ggplot(ciber_fc2, aes(x=FcR, y=Expression, fill=Cell.Type))+
  geom_bar(width = 1, stat = "identity") + labs(title = "Day 7 FcR expression",x="\nFcR", y="") + ylim(0,120000) + theme(axis.text.x = element_text(size= 11,face="bold"))

bp_7


#Original Results ####
og_ciber <- fread("CiberSortInput.txt")
rownames(og_ciber) <- og_ciber$Gene
og_ciber <-as.matrix(og_ciber)
exprs_1 <- og_ciber[,grep("V1_24", colnames(og_ciber))] #Just 12h data

#Group plot ####
list <- list(bp_D0, bp_1, bp_7)
Plot <- wrap_plots(list,ncol = 3,nrow = 1)

Plot + plot_annotation(tag_levels = 'A')  +
  plot_layout(guides = 'collect')

#1600, 750

#No eosinophil/mast cells
bp_1 <-ciber_1fc2 %>% filter(Cell.Type != "Mast cells",Cell.Type != "Eosinophils") %>%
  ggplot(aes(x=FcR, y=Expression, fill=Cell.Type))+
  geom_bar(width = 1, stat = "identity") + labs(title = "Day 1\n", x="\nFcR", y="") + ylim(0,60000) + theme(legend.position = "none", axis.text.x = element_text(size= 11,face="bold")) +scale_fill_manual(values=wes)

bp_1

bp_D0<- ciber_0fc2 %>% filter(Cell.Type != "Mast cells",Cell.Type != "Eosinophils") %>%
  ggplot(aes(x=FcR, y=Expression,fill=Cell.Type))+
  geom_bar(width = 1, stat = "identity") + labs(title = "Baseline\n", x="\nFcR") + ylim(0,60000) + theme(legend.position = "none", axis.text.x = element_text(size= 11,face="bold")) +scale_fill_manual(values=wes)

#+scale_fill_brewer(palette = "Accent")
wes <- c("#46AC8C", "#0B775E","#C6CDF7","#7294d4","#FD6467","#EBCC2A","#E58601", "#B40F20")

bp_D0

bp_7 <-  ciber_fc2 %>% filter(Cell.Type != "Mast cells",Cell.Type != "Eosinophils") %>%
  ggplot(aes(x=FcR, y=Expression, fill=Cell.Type))+
  geom_bar(width = 1, stat = "identity") + labs(title = "Day 7\n",x="\nFcR", y="") + ylim(0,60000) + theme(axis.text.x = element_text(size= 11,face="bold"))+scale_fill_manual(values=wes)

bp_7


list <- list(bp_D0, bp_1, bp_7)
Plot <- wrap_plots(list,ncol = 3,nrow = 1)

Plot + plot_annotation(tag_levels = 'A')  +
  plot_layout(guides = 'collect')
