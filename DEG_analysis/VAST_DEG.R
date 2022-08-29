#Date 25/08/22

#Title - VAST exploratory RNA-seq analysis

#Aims:
#check DEG expression of VAST RNA seq data for before and after vaccination 
#DEG of Fc receptor genes
#Also BAFF/APRILL and CD21
#Eventually set up for cibersort

#Previous script saved in Fc_/RNAseq_setup.R

#Overall aims: ####
#Associate RNAseq expression data with genotype using matrix-eQTL
#Explore FcR expression post-vaccination and challenge
#Associate expression of FcRs with outcome

setwd("~/RNA/Daniel")

#Load packages and data ####
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(tidylog)
library(limma)
library(stringi)
library(janitor)
library(stringr)
library(splines)
library(edgeR)
library(BiocManager)
library(DESeq2)
library(Glimma)
library(ggrepel)

#load(file = "vast_rnaseq.RData")
load(file = "VAST_deg.RData")

#Check data ####

View(x$meta_data) #Same as PData()
View(x$counts) #colnames
View(x$samples)

#Previous steps: (RNAseq_setup.R)

#Normalisation, calculating variance using limma voom
#Checked for batch effects
#Cleaned data of pesky NA samples and matrix dimensions


# DEG: Pre-vac (V1) -- 7 and 28-days post-vaccination (V7, D0) ####

#Voom design matrix ####
str(x$meta_data)
colnames(x$meta_data) <- make_clean_names(colnames(x$meta_data))

#Check factors
x$meta_data$sequence_pool <- as.factor(x$meta_data$sequence_pool)

#Design matrix code
design <- model.matrix(~0 + time_point3 + diagnosis + age_at_do + sex + vaccine_x + sequence_pool,data = x$meta_data)

#Sort and tidy design columns
#514 samples/rows #also extra diagnosis column made
(colSums(design) == 0) #Todo if colsums == 0, print colname, save as vector to filter next
#Remove 0 cols
design <- design[,-c(10:11)]
#Tidy names
colnames(design)<-gsub("vaccine_x","",colnames(design))
colnames(design)<-gsub("time_point3","",colnames(design))
#Non-valid names: D0+12h
colnames(design) <- gsub("\\+", "plus", colnames(design))

#Calculate blocking factor for model and fit ##
corfit <- duplicateCorrelation(v, design=design, block=x$meta_data$participant_id_x)
corfit$consensus #0.33

fit <- lmFit(v, design, block=x$meta_data$participant_id_x, correlation =corfit$consensus)

#Contrast matrix start ##############
cm <- makeContrasts(D0-V0, levels = colnames(design))#TD-D0, #25/08 - D0-V1, V7-V0

#Estimate gene expression values ####
efit <- eBayes(contrasts.fit(lmFit(v,design,block=x$meta_data$participant_id_x,correlation=corfit$consensus), cm))

plotSA(efit) #Norm looks good

#Save all results
D0_V0 <- topTable((efit),num = Inf) #V1_V0
#Check to see that the variance trend has been removed


#Neut genes ####
#BAFF (TNFSF13B), APRIL (TNFSF13), IL-21, IL-17, FcR
#Investigate time course expression if possible 
#Volcano plots and also splines modelling as per FcRs

#Volcano Plot Set up ####

deg = D0_V0 #name of table

p <- ggplot(data=deg, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()
p #Note! skewed volcano were seen before by Daniel

#make genes all lowercase
deg$gene_name <- str_to_lower(deg$gene_name)

#mutate BAFF and APRIL
deg$gene_name <- str_replace(deg$gene_name, "tnfsf13b", "baff")
deg$gene_name <- str_replace(deg$gene_name, "tnfsf13", "april")

#Genes of interest to highlight
gene_list <- c("fcgr3a","fcar","fcgr1a","fcgr1b","fcgr2a","fcgr2b","fcgr2c","april", "baff")

#Make data table with absolute FC values of genes of interest, in this case FCR genes https://www.geeksforgeeks.org/calculate-the-absolute-value-in-r-programming-abs-method/
#For main data mutate a new variable, reg, if FC and P values are above/below a certain threshold

deg <- deg %>%
  mutate(reg =
           case_when(
    deg$logFC >= 0 & deg$adj.P.Val <= 0.05 ~ "Sig Adj. P <0.05",
    deg$logFC <= 0 & deg$adj.P.Val <= 0.05 ~ "Sig Adj. P <0.05",
    deg$logFC >= 0 & deg$P.Value <= 0.05 ~ "Sig P <0.05",
    deg$logFC <= 0 & deg$P.Value <= 0.05 ~ "Sig P <0.05",
    abs(deg$logFC) <= 0 & deg$adj.P.Val >= 0.05 ~ "No Change",
    abs(deg$logFC) <= 0 & deg$adj.P.Val <= 0.05 ~ "No Change",
    abs(deg$logFC) > 0 & deg$adj.P.Val >0.05 ~ "No Change")) %>%
  mutate(reg =
           factor(reg, levels =
                    c("Sig Adj. P <0.05", "Sig P <0.05", "No Change")))

label <- deg[(deg$gene_name %in% gene_list),]

#Plot volcano plot
deg %>% ggplot(aes(x=logFC,y=-log10(P.Value),label=gene_name))+
  geom_point(aes(color = P.Value))+
  labs(title = "28-days Post-Vaccination") +
  theme_minimal() +theme(legend.position = "none") +
  geom_hline(yintercept = 1.4, linetype = 2.5, alpha =0.7) +
  geom_hline(yintercept =2.72, linetype = 2.5, alpha =0.7)+
  geom_label_repel(data=label,size=4,direction="both",nudge_y =1.8,nudge_x =0.1,angle= 70,vjust= 0,segment.size= 0.5,segment.color="#331002",fill="lightgrey")+
  scale_color_gradientn(colours = c("yellow","orange","#ef5a3a","#800000","darkred"),values=c(0,0.0019,0.05,1))


#0 - 0.05 = p-val ns, #work out adj p-value reg-pval sig
#test <- deg %>% filter(adj.P.Val < 0.05) #P =0.0019
#data=label,size=4,direction="both",nudge_y =1.8,nudge_x =0.08,angle= 70,vjust= 0,segment.size= 0.5,segment.color="#331002",fill="lightgrey"
  
degD0V0 <- deg


save.image(file = "VAST_deg.RData")


#Spline regression modelling ####

#Notes####
#Fits model to raw expression values
#but ideally should use estimated expression values from model
#Start with baseline expression values (no modelling)
#Then average exprs values from contrast matrix results tables (also keep logFC)

#V0, V1-V0, V7-V0, D0-V0, Post-challenge time points :)
#Do seperately for TD and nTD participants
#Vi-TCV participants v Vi-PS participants

#Make new data table with these values for genes of interest
#Then plot splines based on average expression values 

#Or can I add covars into splines model as long as I keep them in the ggpot DF
#read up on geom_smooth



exprs <- v$E
gene <- v$genes
target <- v$targets
designV <- v$design

#514 participants

#update colnames by using meta
IDs <- meta$lab_id
Time <- meta$time_point3
test <- paste(IDs,Time,sep = "_") #can use gsub to clean
colnames(exprs) <- test

#could subset data for genes of interest
#rownames of exprs are genes, colnames are participants
#Remake row names as gene names 

test <- cbind(gene,exprs)






Fc_exprs %>% ggplot(aes(x=time, y=fcgr2a, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14) + ylim(10,13.3) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 3)) +
  xlab("Time (days)") + ylab(expression(paste("Fc",gamma,"R2",alpha, " Expression")))
