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

load(file = "vast_rnaseq.RData")

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

#Calculate blocking factor for model ####
corfit <- duplicateCorrelation(v, design=design, block=x$meta_data$participant_id_x)
corfit$consensus #0.33
fit <- lmFit(v, design, block=x$meta_data$participant_id_x, correlation =corfit$consensus)

#Contrast matrix ####
cm <- makeContrasts(V7-V1, levels = colnames(design))#TD-D0, #25/08 - D0-V1, V7-V1


#Estimate gene expression values ####
efit <- eBayes(contrasts.fit(lmFit(v,design,block=x$meta_data$participant_id_x,correlation=corfit$consensus), cm))

plotSA(efit) #Norm looks good

#Save all results
V7_V1 <- topTable((efit),num = Inf)
#Check to see that the variance trend has been removed


#Neut genes ####
#BAFF (TNFSF13B), APRIL (TNFSF13), IL-21, IL-17, FcR
#Investigate time course expression if possible 
#Volcano plots and also splines modelling as per FcRs

#Volcano Plot Set up ####

deg = V7_V1 #name of table

p <- ggplot(data=deg, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()
p #Note! skewed volcano were seen before by Daniel

#make genes all lowercase
deg$gene_name <- str_to_lower(deg$gene_name)

#Genes of interest to highlight
gene_list <- c("fcgr3a","fcar","fcgr1b","fcgr2a","il21ra","tnfsf13", "tnfsf13b","il17ra")


#Make data table with absolute FC values of genes of interest, in this case FCR genes https://www.geeksforgeeks.org/calculate-the-absolute-value-in-r-programming-abs-method/


#For main data mutate a new variable, reg, if FC and P values are above/below a certain threshold
deg <- deg %>%
  mutate(reg = case_when(
    deg$logFC >= 0 & deg$adj.P.Val <= 0.05 ~ "UP",
    deg$logFC <= 0 & deg$adj.P.Val <= 0.05 ~ "DOWN",
    abs(deg$logFC) <= 0 & deg$adj.P.Val >= 0.05 ~ "no_change",
    abs(deg$logFC) <= 0 & deg$adj.P.Val <= 0.05 ~ "no_change",
    abs(deg$logFC) > 0 & deg$adj.P.Val >0.05 ~ "no_change"
  )) %>%
  mutate(reg = factor(reg, levels = c("UP", "no_change","DOWN")))

#Plot volcano plot 
deg %>% ggplot(aes(x=logFC,y=-log10(P.Value), color = reg)) + geom_point()

#Set up labels 
label <- deg[(deg$gene_name %in% gene_list),]

#Plot Volcano ####
deg %>% ggplot(aes(x=logFC,y=-log10(P.Value),label=gene_name))+
  geom_point(aes(color=adj.P.Val))+
  scale_color_gradientn(colours = c("#a5342d","darkred", "orange", "yellow"),values=c(0,0.011,1))+
  theme_minimal()+
    geom_label_repel(data=label,size=4,direction="y",nudge_y =4,nudge_x =-0.15,angle= 60,vjust= 0,segment.size= 0.5,segment.color="black",fill="grey")+
  labs(title = "7-days Post-Vaccination")

??geom_label_repel




               
