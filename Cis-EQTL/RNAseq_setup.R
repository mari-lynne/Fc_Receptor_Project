#RNA seq analysis

#Overall aims: ####
#Associate RNAseq expression data with genotype using matrix-eQTL
#Explore FcR expression post-vaccination and challenge
#Associate expression of FcRs with outcome

#Steps#
# 1) Explore the data
# 2) Format for matrix-eqtl
# 3) Format for DEG analysis 

setwd("~/RNA/Daniel")

#Packages ####
#load packages ####
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

load("Filter2_VAST_STAR_HTSeq_gene_meta_autosomes_mismatch_corrected_demo_minus_rRNA_globins_autosomes_2021-04-27.R")

#load("Filter2_TyGER_combined_RNAseq_runs_STAR_HTSeq_autosome_gene_meta.HLA_genotype_corrected_minus_rRNA_globins_demo_autosomes_2021-05-18.R")

samp <- VAST_autosomes$samples
exprs <- VAST_autosomes$counts
symbol <- VAST_autosomes$gene
pData <- VAST_autosomes$meta_data

#Explore DEG data ####
#Tidy Pheno Data ####
pData <- clean_names(pData)
#Rename days since challenge var
names(pData)[names(pData) == "days_since_challenge"] <- "time"

#Update row names of exprs matrix with gene names
genes <- symbol$gene_name
rownames(exprs) <- genes
#Modify duplicate gene names
rownames(exprs) <- make_clean_names(rownames(exprs))

#Instead of all the transposing lets just subset pData
#Then use row_names from pData to select/subset cols 

#First time point to investigate 0-12h
#Contrast filter (change to suit contrasts)
pData2 <- filter(pData, time == "0"| time == "0.5") #205 rows
pData2 <- pData2 %>% drop_na(participant_id_x)

exprs2 <- exprs[,colnames(exprs) %in% pData2$row_names] #205 cols :)

#Make Design and contrast matrix ####

pData2$time <- as.factor(pData2$time)

design_full <- model.matrix(~0 + time + diagnosis + age_at_do + sex + vaccine_x, data = pData2)

design_full <- design_full[,-c(4:5)]
#203 rows - we have 205 rows in exprs data so they don't match

#If doing a time point DEG/fold change then we need to include participant into the model
#This is done by factoring in the correlation

corfit <- duplicateCorrelation(exprs2,design=design_full,block=pData2$lab_id)

corfit$consensus #Then this inter-subject correlation is input into the linear model fit:
fit <- lmFit(exprs2,design=design_full,block=pData2$lab_id,correlation=corfit$consensus)

#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm <- makeContrasts(time0.5-time0, levels = colnames(design_full))


ebayesfit <- eBayes(contrasts.fit(lmFit(exprs2,design=design_full,block=pData2$lab_idr,correlation=corfit$consensus), cm))


results <- topTable((ebayesfit),num = Inf)
library(ggplot2)
p <- ggplot(data=results, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()

p
#Data doesn't seem to be normalised 
