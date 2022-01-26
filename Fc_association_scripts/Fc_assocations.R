#FcY association testing ####

setwd("~/GWAS/all_enteric/clean_data/post-impute/Plink")
# FcY Range :) -chr 1 --from-bp 1161504430 --to-bp 161729143

#GLM uses logistic regression model
#The 'firth-fallback' modifier requests logistic regression, followed by Firth regression whenever the logistic regression fails to converge. This is now the default.
#1 binary phenotype loaded (173 cases, 139 controls).
#9 covariates loaded from pca_covar.txt.


#No LD pruning ####
#Includes Typhoid/Paratyphoid and Vaccinated individuals

system("./plink2 --pfile enteric_QC --chr 1 --from-bp 161504430 --to-bp 161729143 --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4 --covar-variance-standardize --glm --out enteric_assoc")

gwasResults <- fread("enteric_assoc.Diagnosed.glm.logistic.hybrid")

gwasResults$TEST <- as.factor(gwasResults$TEST)

gwasResults <- filter(gwasResults, TEST == "ADD")

#487 SNPs tested in region

#Rename cols
colnames(gwasResults) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1","FIRTH?","TEST","OBS_CT","OR", "LOG(OR)_SE","Z_STAT","P","ERRCODE")

#Add bonferrnoni column

gwasResults <- mutate(gwasResults, Bonf = P *487) # should be 487 testing 26

#sort by P value
sig <- filter(gwasResults, P <0.05) #73
MTsig <- filter(gwasResults, Bonf <0.05) 

#with the *26 bonferroni correction one snp is significant still
#161722364 chr1:161722364:T:C

#without LD pruning no significant results under multiple testing :(

#6 SNPs with P <0.01
sig <- filter(gwasResults, P <0.001)

#Basic LD pruning and assoc test ####

system("./plink2 --pfile enteric_LD --chr 1 --from-bp 161504430 --to-bp 161729143 --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --out Fc_LD")

#26 Varients

gwasResults <- fread("Fc_LD.Diagnosed.glm.logistic.hybrid")

gwasResults$TEST <- as.factor(gwasResults$TEST)

gwasResults <- filter(gwasResults, TEST == "ADD")

#Rename cols
colnames(gwasResults) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1","FIRTH?","TEST","OBS_CT","OR", "LOG(OR)_SE","Z_STAT","P","ERRCODE")

gwasResults <- mutate(gwasResults, Bonf = P *26)

#sort by P value
sig <- filter(gwasResults, P <0.05) #4
MTsig <- filter(gwasResults, Bonf <0.05) # 0

#without LD pruning no significant results under multiple testing :(

#2 SNPs with P <0.01
sig <- filter(gwasResults, P <0.01)

#Adding adjust flag ####

system("./plink2 --pfile enteric_LD --chr 1 --from-bp 161504430 --to-bp 161729143 --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --adjust --out Fc_LD")

#genomic inflation super high 3.5 I will need to fix this :/
#acc makes sense because there are a lot of sig snps more so than you'd expect this isnt the genomic inflation factor of the whole genome

gwasResults <- fread("Fc_LD.Diagnosed.glm.logistic.hybrid.adjusted")


#Notes ####
#So far these results are randomly pruned for linkage disequilibrium, could miss out on the more significant SNPs see snp_annotating script and haploview for extensions

#For now let's visualize our top snps with no pruning/prioritisation

#SNP Case-Control Plotting ####

#Do top 'protective' SNP and top suscepitible SNP
#Case control analysis 


#Filter original PED file (this has genotype information for SNPs of interest)
#Merge with pheno file
#Plot SNP genotype by cases and controls :)

#Steps

#Update genotype information to RsID ####

#see pos2rsID script
#Use Fcy_rs pgen files



#make Fc Ped files

system("./plink2 --pfile enteric_LD --chr 1 --from-bp 161504430 --to-bp 161729143 --make-bed --out Fcy")

system("./plink --bfile Fcy --recode --out Fcy")

  bim <- fread("Fcy.bim")
  
#try get frequency from --fisher
  #turn frequency into counts mutate 

system("./plink2 --pfile enteric_QC --chr 1 --from-bp 161504430 --to-bp 161729143 --make-pgen --out Fcy") #487 varients

system("./plink2 --pfile enteric_QC --chr 1 --from-bp 161504430 --to-bp 161729143 --make-bed --out Fcy")

system("./plink2 --pfile Fcy --geno-counts")

counts <- fread("plink2.gcount")

system("./plink --bfile Fcy --pheno pheno_plink1.txt --pheno-name Diagnosed --allow-no-sex --assoc counts")

counts <- fread("plink.assoc")

system("./plink --bfile Fcy --pheno pheno_plink1.txt --pheno-name Diagnosed --allow-no-sex --assoc")

freq <- fread("plink.assoc")

#Allele 1 is usually minor
#C_A	Allele 1 count among cases
#C_U	Allele 1 count among controls

sig_counts <- filter(counts, P < 0.01)

#Merge Fc Ped with pheno_update.txt
#Filter for sig SNPs
#Recode poss
#Make stacked bar chart
ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity")

#try use original ped file? #wont be merged
#output from imputation was vcf.gz, merged files then converted to bed > pgen


#convert enteric_QC.bed

system("./plink --bfile enteric_QC --chr 1 --from-bp 161504430 --to-bp 161729143 --make-bed --out Fcy") #--recode --tab 
test <- fread("test.ped")


#BEDMatrix ####
install.packages("remotes")
remotes::install_github("QuantGen/BEDMatrix")
library(BEDMatrix)

m2 <- BEDMatrix("~/GWAS/all_enteric/clean_data/post-impute/Plink/test", n = NULL, p = NULL, simple_names = TRUE)

colnames(m2)
colnames(m[,5:10])
rownames(m2)

m[1:4,400] #what are 0 1 2 encodings

#in the data it is obvious that C is the major allele, and T is the minor allele. So CC is coded as 0, CT is 1,and TT 2.
#I think 0 is homozygous for minor allele, 1 is hetro, 2 is homozygous

#try plotting for SNP chr1:161722364:T:C #have to add _C
#Major allele is T, Ref/Minor is C

vector <- m[, "chr1:161722364:T:C_C"]

df <- as.matrix(m2)
df <- as.data.frame(df)

#need to left join with phenocolumn
pheno <- fread("pheno_update.txt") #IDs are different 


#Hardy option ####

#plink 1.9, check plink 2
--hardy

system("./plink --bfile enteric_QC --chr 1 --from-bp 161504430 --to-bp 161729143 --make-bed --out Fcy")

system("./plink --bfile Fcy --pheno pheno_plink1.txt --pheno-name Diagnosed --allow-no-sex --hardy") 

#227 phenotype values present after --pheno

hardy <- fread("plink.hwe")

hardy <- hardy[,1:6]

#Plink2 version hardy ####

system("./plink2 --pfile Fcy --pheno pheno_update.txt --pheno-name Diagnosed --hardy --keep-if Diagnosed == 1 --out nTD")#139 controls

nTD <- fread("nTD.hardy")

#use keep if to make a diagnosed table, can then compare to undiagnosed table

system("./plink2 --pfile Fcy_rs --pheno pheno_update.txt --pheno-name Diagnosed --hardy --keep-if Diagnosed == 2 --out TD") #173 cases

TD <- fread("TD.hardy")




#SNPtest ####

#Background ####
#Plink2 does not support visualisation of allele types by individual snp/case control very well
#Previous method used ped files but back conversion from PLINK2 to PED is very tricky and also the code used to visualise ped/pheno data wasn't great either

#Intergrate downstream analysis with SNPtest
#So use PLINK for processing of data and genome/region wide summary statistics
#Pass important/significant snps into SNPtest
#Breaks down SNPs by genotype etc.
#requires bgen format

#Convert to BGEN ####
./plink2 --bfile enteric_QC --chr 1 --from-bp 161504430 --to-bp 161729143 --make-bed --out Fcy

#Github code


