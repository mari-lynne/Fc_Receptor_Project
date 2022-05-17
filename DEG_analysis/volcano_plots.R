#Volcano plots 16/05/22 ####

#Data set up
#Data generated in time_course.R script and deg_may 

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

#Filter nas/make time point df's ####
#remove nas for any analysis variables
#Do at this stage so samples from both the pheno and exprs rows with na's are removed
#Also when you can filter time points
pheno_exprs <- pheno_exprs %>% drop_na(time) %>% filter(time >=0)
#Filter minus time points for now
#Re-split 
pData2 <- pheno_exprs[,c(1:14)]
exprs2 <- pheno_exprs[,c(14:10461)] 


#Contrast filter (change to suit contrasts)
T0_12h <- filter(pheno_exprs, time == "0"| time == "0.5")

pData2 <- T0_12h[,1:14]
exprs2 <- T0_12h[,15:10461]
#480 rows remaining

#Reform eset obeject
#The number of rows in phenoData must match the number of columns in assayData
#Row names of phenoData must match column names of the matrix / matrices in assayData.

#Make Time_ID row names in pheno data
rownames(pData2) <- pData2$time_id

#Reform eset obeject
#re transpose
exprs2 <- t(exprs2)
#convert exprssion data to eset
eset2 <- ExpressionSet(assayData = exprs2)
exprs2 <- sapply(exprs2, as.numeric)

#2) Design Matrices ####
table(pData2$time) 
#Set up vars for design matrix ####
#make time as numeric if modelling across time or factor for direct comps
#The design matrix tells the linear model if a particular condition is on, therefore to calculate/model the expression for that value
#So if comparing TD vs nTD the on off switch will be for that (pre-select timepoint)
#If comparing expression between time points then the on off switch is needed for that too
pData2$time <- as.factor(pData2$time)
pData2$study_arm <- as.factor(pData2$study_arm)
#make new vaccine factor (will give more df)
pData2$vax <- ifelse(pData2$study_arm == "T1",0,
                     ifelse(pData2$study_arm == "Placebo",0,
                            ifelse(pData2$study_arm == "M01ZH09", 1,
                                   ifelse(pData2$study_arm == "Ty21a", 1,NA))))
pData2$vax <- as.factor(pData2$vax)

design_full <- model.matrix(~0 + time + diagnosis + age + sex + vax, data = pData2)

#add blocking var for participant
#Then we estimate the correlation between measurements made on the same subject:
corfit <- duplicateCorrelation(eset2,design=design_full,block=pData2$part_number)

corfit$consensus #Then this inter-subject correlation is input into the linear model fit:
fit <- lmFit(eset2,design=design_full,block=pData2$part_number,correlation=corfit$consensus)

#Now we can make any comparisons between the experimental conditions in the usual way, example:
cm <- makeContrasts(time0-time0.5, levels = colnames(design_full))


ebayesfit <- eBayes(contrasts.fit(lmFit(eset2,design=design_full,block=pData2$part_number,correlation=corfit$consensus), cm))


results <- topTable((ebayesfit),num = Inf)


#Volcano Plot ####
p <- ggplot(data=results, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()

p

#Rename deg results
deg <- results

#Genes of interest to highlight
gene_list <- c("fcgr3a","fcar","fcgr1b","fcrla","fcrlb","fcgr2a")


#Make data table with absolute FC values of genes of interest, in this case FCR genes https://www.geeksforgeeks.org/calculate-the-absolute-value-in-r-programming-abs-method/
data= subset(deg, rownames(deg) %in% gene_list)
#grep labels on row.names
deg$label<-row.names(deg) 


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
deg %>% ggplot(data=resultsaes(x=logFC,y=-log10(P.Value))+ geom_point(aes(color=reg))
               
#Volcano Plot ####
ggplot(data=deg, aes(x=logFC, y=-log10(adj.P.Val),label=label)) + geom_point(aes(color=adj.P.Val))+ scale_color_gradientn(colours = c("#a5342d","darkred", "orange", "yellow"), values=c(0,0.011,1)) + theme_minimal()+
geom_label_repel(data= data,size=4,direction="y",nudge_y =4,nudge_x =-0.15,angle= 60,vjust= 0,segment.size= 0.5,segment.color="black",fill="grey") + labs(title = "Baseline-12h Post-Challenge") #a5342d #b5651d

















#Nested analysis make new timeTD factor
Treat <- factor(paste(pData2$diagnosis,pData2$time,sep="."))
design <- model.matrix(~0+Treat)
colnames(design) <- levels(Treat)


design_full <- model.matrix(~time + diagnosis + age + sex + vax, data = pData2)
#time is currently linear, we want factor comparisons also 
#Contrast matricies ####
#Make Contrasts ####
contrast_matrix <- makeContrasts(time0.5,)

?makeContrasts

ebayesfit <- eBayes(contrasts.fit(lmFit(exprs2, exp_design = exp_design), contrast_matrix))
#Error row names of contrasts don't match col names of coefficients

tableexp <- topTable(ebayesfit, number = Inf)






lm <- lmFit(exprs2, design = design_full)
summary(lm)
fit <- eBayes(lm)
topTable(fit)



#Make design matrix
exp_design <- model.matrix(~Time + Array_experiment + age + Sex + StudyArm + Ancestry, data = exp_design)

#Format for Time Course analysis ####
levels(Fc_exprs$Time)







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


#1) Volcano plot 

contrast_matrix <- makeContrasts(Time-0.5, levels = exp_design)

ebayesfit <- eBayes(contrasts.fit(lmFit(exprs2, exp_design = exp_design), contrast_matrix))

tableexp <- topTable(ebayesfit, number = Inf)