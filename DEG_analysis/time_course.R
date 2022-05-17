#Time course analysis of T1/T2 challenge data


#1) Data set up
#2) Annotation
#3) Design matricies
#4) Contrast matricies - Model fitting
#5) Plot top DEG's and genes of interest
#6) Pairwise comparisons

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

#1) Data set up ####
#Import microarray data ####
load("T1_T2_with_T2_samples_rsn.norm.fil.eSet.Ensembl_ID_AVE_collapsed_autosomes_2021-04-28.R",envir = .GlobalEnv)

pData <- pData(T1_T2_autosomes)
eset <- exprs(T1_T2_autosomes)
myPca <- prcomp(t(eset))
exprs <- as.data.frame(eset)
fdata <- fData(T1_T2_autosomes)

#Tidy Pheno Data ####
pData <- clean_names(pData)
#Make new typhoid diagnosis var
pData <- pData %>% mutate(diagnosis = ifelse(is.na(pData$day_of_typhoid_diagnosis_from_challenge), "nTD", "TD"))
#Rename days since challenge var
names(pData)[names(pData) == "days_since_challenge"] <- "time"

#Make new names (issue with matching Ids previously)
pData <- pData %>% mutate(time_id = paste(part_number, time, study_arm, sep = "_"))
#pData <- pData %>% mutate(new_id = paste(part_number, study_arm, sep = "_"))
pData$time_id <- str_replace(pData$time_id, "-", "min")

#Recode variables
pData$time <- as.numeric(pData$time)
pData$diagnosis <- as.factor(pData$diagnosis)


#2) Update gene annotation ####

symbol <- fdata %>% select(Symbol)
#Modify duplicate gene names
symbol$Symbol <- make_clean_names(symbol$Symbol)

#Rows with ENS ID's currently aren't a column variable so we can't update the tables
#Reformat row names to a new column variable, then transpose df so can be re-joined with pheno data

symbol <- tibble::rownames_to_column(symbol, "E_ID")
table <- tibble::rownames_to_column(exprs, "E_ID")
symbol_exprs <- left_join(symbol,table, by = "E_ID") %>% dplyr::select(-E_ID)

#Transpose data frame 
exprs_t <- transpose(symbol_exprs)
names(exprs_t) <- exprs_t[1,] %>% as.vector(mode = "character")
exprs_t <- exprs_t[-c(1),] 

# 2b) Bind tables ####
pheno_exprs <- cbind(pData, exprs_t)
rm(symbol, table, exprs_t)
#symbol_exprs is in the format needed for limma later

#remove nas for any analysis variables
#Do at this stage so samples from both the pheno and exprs rows with na's are removed
pheno_exprs <- pheno_exprs %>% drop_na(time) %>% filter(time >=0)
#Filter minus time points for now 

#Re-split 
#10461 vars - #14 vars of pData
10461-14

pData2 <- pheno_exprs[,c(1:14)]
exprs2 <- pheno_exprs[,c(14:10461)]

#Reform eset obeject
#The number of rows in phenoData must match the number of columns in assayData
#Row names of phenoData must match column names of the matrix / matrices in assayData.

#Make Time_ID row names in pheno data
rownames(pData2) <- pData2$time_id 


#Re Transpose exprs data frame 
exprs2 <- transpose(exprs2)
names(exprs2) <- exprs2[1,] %>% as.vector(mode = "character")
exprs2 <- exprs2[-c(1),] 

#recode to numeric
exprs2 <- sapply(exprs2, as.numeric)


#eset2 <- ExpressionSet(phenoData = phenoData(pData2),cassayData = exprs(exprs_t))

#2) Design Matrices ####
#blocking factor splines (research later I think it's not necessary)
table(pData2$time) #Should I also keep in timepoint3 data and diagnosis

#Set up vars for design matrix ####
#make time as numeric
pData2$time <- as.numeric(pData2$time)
pData2$study_arm <- as.factor(pData2$study_arm)
#make new vaccine factor (will give more df)
pData2$vax <- ifelse(pData2$study_arm == "T1",0,
                          ifelse(pData2$study_arm == "Placebo",0,
                               ifelse(pData2$study_arm == "M01ZH09", 1,
                                      ifelse(pData2$study_arm == "Ty21a", 1,NA))))
pData2$vax <- as.factor(pData2$vax)
#Make design matrix ####
#can be done manually by selecting appropriate vars as a df
#or using model.matrix function
#Check matrix against design on paper 

design_full <- model.matrix(~diagnosis + time + age + sex + vax, data = pData2)


#Splines ####
library(splines)

#We also need to add splines points to our design matrix so lm function can use these to estimate the model
#Make a separate splines matirx or df before combining this with our original design

design_splines <- ns(pData2$time, df = 5)

#Full matrix ####
#design_splines_full <- model.matrix(~design_full + design_splines)
#Just bind matricies as they're compatible
design_splines_full <- cbind(design_full, design_splines)

#Model per diagnosis across time using splines ####
library(limma)
lm <- lmFit(exprs2, design = design_splines)

summary(lm)#looks weird imo

fit <- eBayes(lm)
top <- topTable(fit)

top










#Test effect of just time (nTD and TD)
# Create a design matrix to test the effect of time only, whatever the treatment
design_limma_time <- model.matrix(~ design_full$diagnosis + design_full)



#Test effect of time just nTD

#Effect of time nTD

#Create a design matrix to test the effect of time and diagnosis


#Measure the effect of diagnosis
X <- ns(targets$time, df=5) #3 splines around each time point value
#Then fit separate curves for the control and treatment groups:
Group <- factor(targets$diagnosis)
design <- model.matrix(~Group*X)


design_splines_full <-  ns(design_full$Time, df = 3)
design_splines_RAS  <-  ns(design_RAS$Time, df = 5)

# Create a design matrix to test the effect of time on RAS timecourse only
design_limma_RAS <- model.matrix(~ design_splines_RAS) #They just filtered for RAS

# Create a design matrix to test the effect of time only, whatever the treatment
design_limma_time <- model.matrix(~ design_full$Treatment + design_splines_full)

# Create a design matrix to test the effect of time and treatment
design_limma_full <- model.matrix(~ 0 + design_full$Batch + design_splines_full : design_full$Treatment)
colnames(design_limma_full) <- gsub(":", ".", colnames(design_limma_full))
colnames(design_limma_full) <- gsub("\\$", ".", colnames(design_limma_full))





#add participant blocking factor
#eset2
emat <- as.matrix(exprs2)
eset2 <- ExpressionSet(assayData = emat)

corfit <- duplicateCorrelation(eset2,design,block=targets$new_id)
corfit$consensus

#Model expression accounting for participant ID
fit <- lmFit(exprs2,design,block=targets$new_id,correlation=corfit$consensus)
#fit <- lmFit(exprs2, design)
fit <- eBayes(fit)
top <- topTable(fit)








#Extra code #####


# 2b) Update gene names/gene expression ####
symbol <- as.data.frame(fData(T1_T2_autosomes)) %>% select(Symbol)
#Modify duplicate names
symbol$Symbol <- make_clean_names(symbol$Symbol)
#Rows with ENS IDs currently aren't a column variable so we can't update, reformat :)
symbol <- tibble::rownames_to_column(symbol, "E_ID")
tableexp <- tibble::rownames_to_column(exprs, "E_ID")
symbol_exprs <- left_join(symbol,tableexp, by = "E_ID") %>% dplyr::select(-E_ID)
#update colnames
exprs_t <- transpose(symbol_exprs)
gnames <- exprs_t[1,] %>% as.vector(mode = "character")
names(exprs_t) <- gnames
exprs_t <- exprs_t[-c(1),]
exprs_t <-sapply(exprs_t, as.numeric)
                                     

#Error in lmFit(exprs_t, design) : 
#row dimension of design doesn't match column dimension of data object
#nrs in design = 572
#ncols exprs = 575??
#y = exprs data
#needed to check for nas in design factors as didn't all make groups     




cor <- duplicateCorrelation(exprs2, design, block=new_id)

cor$consensus.correlation

#get expression data
fit <- lmFit(exprs2, design)
fit <- eBayes(fit)
topTable(fit)

#update enames
enames = colnames(eset)
pnames = pData$time_id

check <- cbind(enames, pnames)
colnames(eset) <- pnames

#remove nas for any analysis variables
pData2 <- select(pData,new_id, time_id, time, diagnosis) %>% drop_na()

exprs2 <- as.data.frame(eset)
exprs2 <- exprs2[,(colnames(eset) %in% pData2$time_id)]

?ExpressionSet


exprs2 <- as.data.frame(eset)
exprs2 <- exprs2[,(colnames(eset) %in% pData2$time_id)]