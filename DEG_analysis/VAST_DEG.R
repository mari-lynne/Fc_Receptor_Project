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

plotSA(efit)

#Save all results
V7_V1 <- topTable((efit),num = Inf)
#Check to see that the variance trend has been removed


#Neut genes ####

#Volcano plots and also splines modelling as per FcRs

#BAFF, APRIL, IL-21, IL-17, FcR
#Investigate time course expression if possible 
