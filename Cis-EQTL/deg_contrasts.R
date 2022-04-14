#Set up T0_12h gene expression data set
library(tidyr)
library(tidylog)
library(dplyr)
library(data.table)

setwd("~/RNA/Daniel")
load("T1_T2_with_T2_samples_rsn.norm.fil.eSet.Ensembl_ID_AVE_collapsed_autosomes_2021-04-28.R",envir = .GlobalEnv)

#Set up Pheno Data ####
pData <- pData(T1_T2_autosomes)
#make new typhoid diagnosis variable
pData <- pData %>% mutate(Diagnosis = ifelse(is.na(pData$Day.of.typhoid.diagnosis..from.challenge), 0, 1))
names(pData)[names(pData) == "Days.since.challenge"] <- "Time"

T0_12h <- T1_T2_autosomes[pData$Time == "0"| pData$Time == "0.5"]

#mSet up Expression Data ###
exprs <- as.data.frame(exprs(T1_T2_autosomes))

#update gene names ####
names <- as.data.frame(fData(T1_T2_autosomes))
symbol <- select(names, Symbol)                    
#Rows with ENS IDs currently aren't a column variable so we can't update, reformat :)
symbol <- tibble::rownames_to_column(symbol, "E_ID")
tableexp <- tibble::rownames_to_column(exprs, "E_ID")
symbol_exprs <- left_join(symbol,tableexp, by = "E_ID") %>% select(-E_ID)

#transpose data frame
df_t <- transpose(symbol_exprs)
gnames <- df_t[1,] %>% as.vector(mode = "character") #update colnames
names(df_t) <- gnames
df_t <- df_t[-c(1),] 

#Bind tables ####
pheno_exprs <- cbind(pData, df_t)

#Filter
T0_12h <- filter(pheno_exprs, Time == "0"| Time == "0.5")

#Prepare data for DEG analysis ####
pData2 <- T0_12h[,1:12]
exprs2 <- T0_12h[,13:10459]

#exprs 2 needs to be converted into a matrix 
#re transpose
exprs2 <- t(exprs2)
#convert exprssion data to eset
eset2 <- ExpressionSet(assayData = exprs2)

#table expsign Matrix ####

pData2$Time<- as.factor(pData2$Time)
pData2$ArrayArray <- as.factor(pData2$Array_experiment)
pData2$Ancestry <- as.factor(pData2$Ancestry)
pData2$StudyArm <- as.factor(pData2$StudyArm)

tableexpsign_12h <- select(pData2, Time, Array_experiment, age, Sex, StudyArm)

#make matrix
tableexpsign <- motableexpl.matrix(~Time + Array_experiment + age + Sex + StudyArm, data = tableexpsign_12h)

#Make Contrasts ####
contrast_matrix <- makeContrasts(Time0.5, levels = tableexpsign)
ebayesfit <- eBayes(contrasts.fit(lmFit(exprs2, tableexpsign = tableexpsign), contrast_matrix))
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


# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
tableexp$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
tableexp$diffexpressed[tableexp$logFC > 0.6 & tableexp$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
tableexp$diffexpressed[tableexp$logFC < -0.6 & tableexp$adj.P.Val < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=tableexp, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + geom_point() + theme_minimal()
p

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
p2



#Results
library(stringr)
library(dplyr)
FCGR <- results %>% filter(str_detect(Symbol,'FCR|FCG'))

#= 12 genes :)
#At 12h no diff expression tho
#Try 1 day and also 24hours 

#need to get expression by participant by gene

#Filter pheno_exprs table for 12 hour data

twelve_hr <- pheno_exprs %>% filter(Time == "0.5") 

