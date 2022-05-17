#Date 09/05/22


#Aims ####

#Set up microarray data from T1T2 typhoid challenge studies
#Create desigin and contrast matricies of interest

#Import data ####

#load packages ####
library(dplyr)
library(data.table)
library(tidyr)
library(tidylog)
library(limma)
library(stringi)

#Functions ####

#Data cleaning 
rownames2col <- function(df, colname) {
  output <- cbind(row.names(df), df)
  colnames(output)[1] <- colname
  return(output)
}
col2rownames <- function(df, colname, removecol=FALSE){
  row.names(df) <- df[,colname]
  if(removecol){df[,colname] <- NULL}
  return(df)
}

'%!in%' <- function(x,y)!('%in%'(x,y))


#Import microarray data ####
setwd("~/RNA/Daniel")
load("T1_T2_with_T2_samples_rsn.norm.fil.eSet.Ensembl_ID_AVE_collapsed_autosomes_2021-04-28.R",envir = .GlobalEnv)

#Summarise data ####
#Use biobase S4 summary functions to access data from eset object

pd <- pData(T1_T2_autosomes)
eset <- exprs(T1_T2_autosomes)
myPca <- prcomp(t(eset))
exprs <- as.data.frame(eset)


#Tidy Pheno Data ####
pd <- clean_names(pd)
#Make new typhoid diagnosis var
pd <- pd %>% mutate(diagnosis = ifelse(is.na(pd$day_of_typhoid_diagnosis_from_challenge), 0, 1))
#rename days since challenge var
names(pd)[names(pd) == "days_since_challenge"] <- "time"

#Update gene names ####

#Update eset object by joining with names data from phenoData
#Names contains both the ESN ID and gene symbol annotation (which we want)
#First need to transform ESN row names to columns for both objects
#Then we can left join the expression data to the names data
#After joining we can rename the first col to contain new gene names before formating back into eset data

symbol <- as.data.frame(fData(T1_T2_autosomes)) %>% select(Symbol)
#Symbol has useful gene names info    
#Update rownames to column for ENSIDs
symbol <- tibble::rownames_to_column(symbol, "E_ID")
table_exp <- tibble::rownames_to_column(exprs, "E_ID")
symbol_exprs <- left_join(symbol,table_exp, by = "E_ID") %>% select(-E_ID)

#Modify duplicate names
symbol_exprs$Symbol <- make_clean_names(symbol_exprs$Symbol)
#colnames back to rownames
rownames(symbol_exprs) <- symbol_exprs$Symbol 
symbol_exprs <- symbol_exprs[,-1]
#Transform expression data back into matrix format for later
eset <- as.matrix(symbol_exprs)

##########################################

#Set up Design and Contrast matrices ###

#As there are multiple samples per person we end up with an incorrect number of levels/comparisons unless we make a filtered dataset (both pData and exprs), and run pairwise comparisons from there eg. pData_0_12h (1)

#We can also model expression across time using linear regression/splines (2)

#1) Data subsetting for pairwise comparisons ####

#Filter
pd0_12h <- filter(pd, time == "0"| time == "0.5") #185
exprs0_12 <- exprs[,names(exprs) %in% pd0_12h$sample_id] #91 vars

#pd will have duplicate entries for each sample_id
#update exprs to have part number IDs instead?


#Prepare data for DEG analysis ####




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



#Need to filter pData and eset/exprs data
pd0_12h <- filter(pd, time == "0"|time == "0.5")

#Select exprs colnames by the study_ID of filtered table #[rows,cols]
exprs0_12 <- exprs[,names(exprs) %in% pd0_12h$sample_id] #91 vars = 90 participants

#Design matrix 12h ####
design <- pd0_12h %>% select(time, age, sex, study_arm)
design <- model.matrix(~time + age + sex + study_arm, data = design)
#Clean levels #colnames(design) <- str_replace(colnames(design), "time","")

#Make Contrasts ####
#contrast_matrix <- makeContrasts(levels=colnames(design)

contrast_matrix <- makeContrasts(time0.5, levels = design)

#Fit model to estimate average gene expression
                                 ebayesfit <- eBayes(contrasts.fit(lmFit(exprs, design = mat), contrast_matrix))
                                 
                                 
                                 #Error row names of contrasts don't match col names of coefficients





# General Design Matrix ####
pd$time<- as.character(pd$time)
#Clean time var to remove - as gives error in contrast levels, non-valid names
pd$time <- str_replace_all(pd$time, "-","minus") #Also try as numeric later
pd$array_experiment <- as.factor(pd$array_experiment)
pd$ancestry <- as.factor(pd$ancestry)
pd$study_arm <- as.factor(pd$study_arm)
pd$part_number <- as.character(pd$part_number)

design <- pd %>% select(time, age, sex, study_arm)
#Make Design matrix ####
design <- model.matrix(~time + age + sex + study_arm, data = design)
#Clean levels
colnames(design) <- str_replace(colnames(design), "time","")


#Make Contrasts ####

contrast_matrix <- makeContrasts(levels=colnames(design)

contrast_matrix <- makeContrasts(time0.5, levels = mat)
ebayesfit <- eBayes(contrasts.fit(lmFit(exprs, design = mat), contrast_matrix))


#Error row names of contrasts don't match col names of coefficients

 topTable(ebayesfit, number = Inf)

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

twelve_hr <- pheno_exprs %>% filter(time == "0.5")