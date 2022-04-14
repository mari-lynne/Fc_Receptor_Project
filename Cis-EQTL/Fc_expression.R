library(limma)


#Plan
#Extract data for Fc Receptors
#Model their expression over time


#Get gene expression data ####
setwd("~/GWAS_22/Fc_receptor/data/")
table <- fread("gene_exprs_T1T2.txt")#from array set up script



Fc_exprs <- table[,grep("^FCR|FCG|FCAR", colnames(table), value = TRUE,)]
Fc_exprs <- table %>%
  dplyr::select(Fc_exprs)
Fc_exprs <- cbind(pData, Fc_exprs)

#table exp_sign Matrix ####
Fc_exprs$Time<- as.factor(Fc_exprs$Time)
Fc_exprs$ArrayArray <- as.factor(Fc_exprs$Array_experiment)
Fc_exprs$Ancestry <- as.factor(Fc_exprs$Ancestry)
Fc_exprs$StudyArm <- as.factor(Fc_exprs$StudyArm)
Fc_exprs$Sex <- as.factor(Fc_exprs$Sex)

exp_design <- select(Fc_exprs, Time, Array_experiment, age, Sex, StudyArm, Ancestry)

#Make design matrix
exp_design <- model.matrix(~Time + Array_experiment + age + Sex + StudyArm + Ancestry, data = exp_design)


#Format for Time Course analysis ####
levels(Fc_exprs$Time)

#For Fc expression Data run a PCA
#Colour by time point

#Can do either a per time point analysis baseline to 12, 24h and TD
#Or a global analysis 
#For vaccine samples could also do a pre-challenge (post-vac analysis), and a post-challenge analysis

#Baseline samples would be better to compare pre-vaccination T2
#Then filter other - time points

#T2 baseline
test <- filter(Fc_exprs, Time == "-32") #10 people = V0 timepoint in TY21a
test <- filter(Fc_exprs, Time == "-30") #20 ppl M01ZH09 and Placebo V0
test <- filter(Fc_exprs, Timepoint3 == "V0")
#Use V0 samples
#T1 baseline = I think Day 0
test <- filter(Fc_exprs, Timepoint3 == "D0" & StudyArm == "T1")

#Recode to be labeled baseline Timepoint3
#But also should I recode the numerical time to = 0 once other time points are removed
#That way can plot global analysis
#Remove other - Time points for now for Day 0 
Fc_exprs$Time <- recode(Fc_exprs$Time, `-32` = "Baseline", `-30` = "Baseline")
Fc_exprs$Time[Fc_exprs$StudyArm == "T1" & Fc_exprs$Time == "0"] <- "Baseline"

#remove other pre-vac time-points
'%!in%' <- function(x,y)!('%in%'(x,y))
Fc_exprs <- Fc_exprs %>% filter(Time %!in% c("-14","-21","-26","-28"))

#Get Packages
BiocManager::install(version="devel")
BiocManager::valid()
BiocManager::install("timecoursedata")
BiocManager::install("moanin")
# From Github
library(moanin)
library(timecoursedata)
# From CRAN
library(NMF)
library(ggfortify)
# From Bioconductor
library(topGO)
library(biomaRt)
library(KEGGprofile)





Fc_exprs$Time <- as.factor(Fc_exprs$Time)
levels(Fc_exprs$Time)

Fc_exprs$Time <- as.character(Fc_exprs$Time)


#Make Contrasts ####
contrast_matrix <- makeContrasts(Time-0.5, levels = exp_design)

ebayesfit <- eBayes(contrasts.fit(lmFit(exprs2, exp_design = exp_design), contrast_matrix))
#Error row names of contrasts don't match col names of coefficients

tableexp <- topTable(ebayesfit, number = Inf)

#This is fold change in expression between time points yay

#Volcano Plot ####
library(ggplot2)
p <- ggplot(data=tableexp, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()

# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

p2
