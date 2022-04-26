# eQTL analysis ###

### Aims ###
# Prepare data for eQTL analysis
# Associate sig fc snps with fc gene expression and outcome
# Visualise results
# Do for microarray and rna-seq data sets separately (before working on data harmonisation)

library(data.table)
library(dplyr)
library(stringr)

### Microarray data ###

setwd("~/RNA/Daniel")

load("T1_T2_with_T2_samples_rsn.norm.fil.eSet.Ensembl_ID_AVE_collapsed_autosomes_2021-04-28.R",envir = .GlobalEnv)

### Check data ###
#https://www.rdocumentation.org/packages/Biobase/versions/2.32.0/topics/ExpressionSet
#https://kasperdanielhansen.github.io/genbioconductor/html/ExpressionSet.html

featureNames(T1_T2_autosomes)[1:10]
sampleNames(T1_T2_autosomes)[1:10]

tail(pData(T1_T2_autosomes)) #pheno data with sample IDs and timepoints
colnames(pData(T1_T2_autosomes))

#Subsetting has two dimensions; the first dimension is genes and the second is samples. It keeps track of how the expression measures are matched to the pheno data

T1_T2_autosomes[,1:10]#10447 features
T1_T2_autosomes[1:10,] #575 features (as there are multiple time points)


gene_names <- featureNames(T1_T2_autosomes)

gene_names <- as.data.frame(gene_names) #ENS transcripts 

#https://www.bioconductor.org/packages/devel/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf

#Other names in featureData > Data: Transcript, ILMN_Gene, RefSeq_ID

featureNames(T1_T2_autosomes)
fvarMetadata(T1_T2_autosomes)

#Make metadata/ pheno data frame ####
meta_names <- colnames(pData(T1_T2_autosomes))
#check classes
sapply(pData(T1_T2_autosomes), class)

#must use labelDescription for biomart package to work
metadata <- data.frame(labelDescription= meta_names, row.names= meta_names) #could rename cols in row.names if I wanted

#AnnotatedDataFrame that conveniently stores and manipulates the phenotypic data and its metadata in a coordinated fashion. Create and view an AnnotatedDataFrame instance with:
phenoData <- new("AnnotatedDataFrame", data=pData(T1_T2_autosomes), varMetadata=metadata)

head(pData(phenoData))

pheno <- data.frame(pData(phenoData))

#Make gene expression df ####

#Ultimately need to merge the phenodata table with the express data table
#Then add in the feature names 
#Then filter

#Extract expression data frame
#In assayData? - extract using exprs function

head(exprs(T1_T2_autosomes))
#Samples are columns, transcripts are rows
exprs <- data.frame(fData(T1_T2_autosomes), exprs(T1_T2_autosomes))

  exprs_names <- names(exprs)
exprs <- select(exprs, -Probe_Id, -Species, -Source, -Search_Key, -Probe_Sequence, -Probe_Coordinates, -Cytoband, -Accession, -Probe_Start, -Array_Address_Id, -Unigene_ID)

#rename cols remove X
#think of it as keeping the letters after the pattern

ID_cols <- str_extract(exprs_names, regex("(?<=X)[\\w]+"))
ID_cols <- na.omit(ID_cols)
first_cols <- exprs_names[1:30]
exprs_names <- c(first_cols, ID_cols)
names(exprs) <- exprs_names

#Transpose so can join to pheno_data
exprs2 <- exprs[,30:594]
exprs2 <- t(exprs2)
names(exprs2) <- exprs2[1,]
exprs2 <- as.data.frame(exprs2)
exprs2 <- exprs2[2:605,]

#make rows column
exprs2 <- tibble::rownames_to_column(exprs2, "SampleID")
pheno_express <- left_join(pheno, exprs2, by = "SampleID")

#Remove time point from ID
T1 <- filter(pheno_express, StudyArm == "T1")

names <- (T1$SampleID)

#str_extract(names, regex("(?<=_)[\\w]+"))
#need to do inverse of this 
T1_IDs <- str_extract(names, regex("[\\d]+(?=_)")) #positive look ahead

T1$SampleID <- T1_IDs

# T2 sorting 
T2 <- filter(pheno_express, StudyArm != "T1")

names <- (T2$SampleID)
#need to search for 3 or 4 consecutive numbers
#P{m,n}	Between m and n instances of P
T2_IDs <- str_extract(names, regex("([\\d]{2,4})"))
T2$SampleID <- T2_IDs #double check T2 Ids with part number
pheno_express <- rbind(T1, T2)

#Use Part number instead #####
pheno_express <- select(pheno_express, - SampleID)
pheno_express <- rename(pheno_express, Sample_ID = Part_number)

#write.table(pheno_express, file = "~/GWAS_22/Fc_receptor/data/pheno_express_T1T2.txt", row.names = F, quote = F, sep = "\t")

table <- fread("~/GWAS_22/Fc_receptor/data/pheno_express_T1T2.txt")

#Format for cybersort/matrixEQTL ####

#Make contrast matrix DEG

library(limma)
library(edgeR)
table

#Make targets file for design matrix

#This is essentially the metadata 





#See the pattern in FcGR expression over time
