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
ggplot(data=results, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal

#Normalise read counts ####
#Weird volcano plots/results
#Data is not normalised/filtered

x <- VAST_autosomes
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
?cpm

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)#median/mean log2 counts per million

#Remove lowly expressed genes ####
table(rowSums(x$counts==0)==521)
#samples with row sum = 0 across all samples (no samples) 
#13609 samples

#Filter express genes
keep.exprs <- filterByExpr(x)
#w All samples appear to belong to the same group
x <- x[keep.exprs, keep.lib.sizes=FALSE]
?filterByExpr
dim(x)
#Now is 12969

#Plot expression levels ####
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x) #redo using colour ramp
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
#lcpm is unormalised, x has been updated/normalised
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
#Not a big change


#Normalise filtered data ####
x <- calcNormFactors(x, method = "TMM")
length(x$samples$norm.factors) #Also check this is na removed later #currently length = 521

#boxplots
#Define colours
group = as.factor(x$meta_data$Sequence_Pool)
col.group<-group
levels(col.group)<-brewer.pal(nlevels(col.group),"Accent")
col.group<-as.character(col.group)
#raw
par(mfrow=c(1,2))
lcpm <- cpm(VAST_autosomes, log=TRUE)
boxplot(lcpm, las=2, col=col.group, main="")
title(main="A) Raw data",ylab="Log-cpm")
#Norm
lcpm <- cpm(x, log=TRUE)
boxplot(lcpm, las=2, col=col.group, main="")
title(main="B) Normalised data",ylab="Log-cpm")

# MDS ####

#Visually examine the factors to include in you linear modelling
#These could include technical factors such as sequencing lane/batch
#And also time points
#Ideally you would want to see differences in time points where we expect to see DEG's, and less variation between technical factors
#If there is no clustering then that factor is less necessary for the model
#you can also test for two factors together using interaction terms
#group = interaction(x$meta_data$Sequence_Pool,x$meta_data$TimePoint3)

#Sequencing lane 
group = as.factor(x$meta_data$Sequence_Pool)

lcpm<-cpm(x,log=TRUE)
par(mfrow=c(1,2))
col.group<-group
levels(col.group)<-brewer.pal(nlevels(col.group),"Accent")
col.group<-as.character(col.group)
plotMDS(lcpm,labels=group,col=col.group)
title(main= "A) MDS - Sequencing Groups")
#plotMDS(lcpm,labels=group,col = as.numeric(group))
#Quite a strong batch effect
#Will need to account for this technical variation as a covar in any models
#Try pre-normalise for it though

#Group of interest
group = as.factor(x$meta_data$TimePoint3)
#Get log count data
lcpm<-cpm(x,log=TRUE)
col.group<-group
#Set colours
palette_Dark2 <- colorRampPalette(brewer.pal(8, "Set2"))
levels(col.group)<-palette_Dark2(length(unique(x$meta_data$TimePoint3)))
col.group<-as.character(col.group)
#Plot
plotMDS(lcpm,labels=group,col=col.group)
title(main= "B) MDS - Time Points")

#try a pre-norm for comparison
#Sequencing lane 
group = as.factor(VAST_autosomes$meta_data$Sequence_Pool)

lcpm<-cpm(VAST_autosomes,log=TRUE)
par(mfrow=c(1,2))
col.group<-group
levels(col.group)<-brewer.pal(nlevels(col.group),"Accent")
col.group<-as.character(col.group)
plotMDS(lcpm,labels=group,col=col.group)
title(main= "A) MDS - Sequencing Groups")


#Batch correction ####
#Strong batch effects not removed by normalisation
#Use combat to account for this
#or limma
lcpm<-cpm(x,log=TRUE)
lcpm <- limma::removeBatchEffect(lcpm, x$meta_data$Sequence_Pool)

group = as.factor(x$meta_data$Sequence_Pool)
par(mfrow=c(1,2))
col.group<-group
levels(col.group)<-brewer.pal(nlevels(col.group),"Accent")
col.group<-as.character(col.group)
plotMDS(lcpm,labels=group,col=col.group)
title(main= "C) MDS - Sequencing Groups")


#Time point redo
group = as.factor(x$meta_data$TimePoint3)
col.group<-group
#Set colours
palette_Dark2 <- colorRampPalette(brewer.pal(8, "Set2"))
levels(col.group)<-palette_Dark2(length(unique(x$meta_data$TimePoint3)))
col.group<-as.character(col.group)
#Plot
plotMDS(lcpm,labels=group,col=col.group)
title(main= "D) MDS - Time Points")

#Cant batch correct on raw values as gives negative counts
#Therefore just have to include sequencing group in the final model
#https://support.bioconductor.org/p/76099/


#Limma voom ####
#Converts counts data to normalised expression values, accounts for vairanc #Also takes in design matrix

#Plots mean count size (x) against variance (y)
#Lower counts tend to have lower variance
#Therefore voom normalises - flattens this trend
#Variances are calculated by fitting a linear model to the data provided by the design matrix

par(mfrow=c(1,1))
v <- voom(x, plot=TRUE)
#Error in voom(x, plot = TRUE) : Negative counts not allowed
#Batch correction results in negative numbers 

#Voom design matrix ####
colnames(x$meta_data) <- make_clean_names(colnames(x$meta_data))
#Include technical factors such as batch etc.
x$meta_data$sequence_pool <- as.factor(x$meta_data$sequence_pool)
design <- model.matrix(~0 + time_point3 + diagnosis + age_at_do + sex + vaccine_x + sequence_pool,data = x$meta_data)

#514 samples/rows #also extra diagnosis column made
(colSums(design) == 0)
#Remove 0 cols
design <- design[,-c(10:11)]

colnames(design)<-gsub("vaccine_x","",colnames(design))
colnames(design)<-gsub("time_point3","",colnames(design))
#Time is currently a factor

#NA sorting ####

#Unequal dims of design matrix and exprs matrix due to na vars
#Solution is to find these na's, filter them, then use the filtered design matrix to subset the orignal expression data so they're equal.
#Currently 514 rows in design, 521 samples in original set

View(x$meta_data) #Same as PData start from here
View(x$counts) #colnames
View(x$samples) #rownames match up

x$meta_data <- x$meta_data %>% drop_na(participant_id_x) #514 rows remaining
#Filter remaining dataframes to keep equal matrix dimensions
x$counts <- x$counts[,colnames(x$counts) %in% x$meta_data$row_names] 
x$samples <- x$samples[rownames(x$samples) %in% x$meta_data$row_names,]

length(x$samples$norm.factors)
#Check DEGlist - all df's now have 514 dims, Redo design matrix

#New voom ####
#Redo Voom with design matrix
#This way normalisation is taking into account other factors 
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

v <- voom(x,design,plot=TRUE)
v


#New Voom 2 (blocking factor) ####
#Add in participant blocking factor
#Then we estimate the correlation between measurements made on the same subject:
corfit <- duplicateCorrelation(v, design=design, block=x$meta_data$participant_id_x)
corfit$consensus #0.33
fit <- lmFit(v, design, block=x$meta_data$participant_id_x, correlation =corfit$consensus)
#Now we can make any comparisons between the experimental conditions in the usual way
cm <- makeContrasts(TD-D0, levels = colnames(design))
#Error in makeContrasts: The levels must by syntactically valid names in R, see help(make.names).  Non-valid names: D0+12h
colnames(design) <- gsub("\\+", "plus", colnames(design))

#Estimate gene expression values 
efit <- eBayes(contrasts.fit(lmFit(v,design,block=x$meta_data$participant_id_x,correlation=corfit$consensus), cm))

#Save all results
results <- topTable((efit),num = Inf)
#Check to see that the variance trend has been removed
plotSA(efit)

save.image(file = "vast_rnaseq.RData")

#SA plot has a lot of low FC's
#Try use  voomWithQualityWeights() - accounts for sample heterogeneity
#Possibly fliter more lowly expressed genes

#Volcano Plot ####
library(ggplot2)
ggplot(data=results, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_minimal()

#TD vs Day0 (all vaccines mixed in)
#Rename deg results
deg <- results

#Genes of interest to highlight
gene_list <- c("FCGR3A","FCGR3C","FCAR","FCGR1B","FCRLA","FCRLB","FCGR2A", "FCGR2B", "FCGR2C", "FCGR1A")


#Make data table with absolute FC values of genes of interest, in this case FCR genes https://www.geeksforgeeks.org/calculate-the-absolute-value-in-r-programming-abs-method/
data= subset(deg, deg$gene_name %in% gene_list)
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
deg %>% ggplot(data=results(aes(x=logFC,y=-log10(P.Value))+ geom_point(aes(color=reg))
                                                                      
#Volcano Plot ####
 ggplot(data=deg, aes(x=logFC, y=-log10(adj.P.Val),label=gene_name)) + geom_point(aes(color=adj.P.Val))+ scale_color_gradientn(colours = c("#a5342d","darkred", "orange", "yellow"), values=c(0,0.011,1)) + theme_minimal()+geom_label_repel(data= data,size=4,direction="y",nudge_y =4,nudge_x =-0.15,angle= 60,vjust= 0,segment.size= 0.5,segment.color="black",fill="grey") + labs(title = "Baseline - Day of Diagnosis")
#a5342d #b5651d

library(ggrepel)




#Voom linear time ####



