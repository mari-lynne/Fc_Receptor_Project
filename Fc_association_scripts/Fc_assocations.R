#FcY association testing ####

library("data.table")
library("dplyr")
library("ggplot2")

setwd("~/GWAS/all_enteric/clean_data/post-impute/Plink")
# FcY Range :) -chr 1 --from-bp 1161504430 --to-bp 161729143
#RS IDs updated using pos2RsID use Fcy_rs files

#GLM uses logistic regression model.
#1 binary phenotype loaded (173 cases, 139 controls).
#9 covariates loaded from pca_covar.txt.

#No LD pruning ####
#Includes Typhoid/Paratyphoid and Vaccinated individuals
#Need to formally assess covariate weighting still

system("./plink2 --pfile Fcy_rs --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4 --covar-variance-standardize --glm --out Fcy")

gwasResults <- fread("Fcy.Diagnosed.glm.logistic.hybrid")
gwasResults$TEST <- as.factor(gwasResults$TEST)
gwasResults <- filter(gwasResults, TEST == "ADD")
#487 SNPs tested in region

#Rename cols
colnames(gwasResults) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1","FIRTH?","TEST","OBS_CT","OR", "LOG(OR)_SE","Z_STAT","P","ERRCODE")
#Add bonferrnoni column
gwasResults <- mutate(gwasResults, Bonf = P *487)

#sort by P value
sig <- filter(gwasResults, P <0.05) #73
MTsig <- filter(gwasResults, Bonf <0.05) 

# None currently meet MT, however this is without pruning so not accurate filtering. With the *26 post-LD bonferroni correction, one snp is significant still #161722364 chr1:161722364:T:C

#6 SNPs with P <0.01
sig <- filter(gwasResults, P <0.01)

#Basic LD pruning and assoc test ####

#Prune Fcy_Rs
system("./plink2 --pfile Fcy_rs --indep-pairwise 50 5 0.2 --out Fcyrs_LD")

system("./plink2 --pfile Fcy_rs --extract Fcyrs_LD.prune.in --make-pgen --out Fcy_rs_LD")
#31 varients
#LD pruning does not prioritise snps of significance, see snp_annotating.R and association_testing.R for attempts to overcome this

system("./plink2 --pfile Fcy_rs_LD --chr 1 --from-bp 161504430 --to-bp 161729143 --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --out Fc_LD")

gwasResults <- fread("Fc_LD.Diagnosed.glm.logistic.hybrid")
gwasResults$TEST <- as.factor(gwasResults$TEST)
gwasResults <- filter(gwasResults, TEST == "ADD")

#Rename cols
colnames(gwasResults) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1","FIRTH?","TEST","OBS_CT","OR", "LOG(OR)_SE","Z_STAT","P","ERRCODE")

gwasResults <- mutate(gwasResults, Bonf = P *31)

#sort by P value
sig <- filter(gwasResults, P <0.05) #3 (previous LD had 4)
MTsig <- filter(gwasResults, Bonf <0.05)

#2 SNPs with P <0.01
sig <- filter(gwasResults, P <0.01)

#even compared to a different LD run (used whole genome), sig snps have been pruned #come back to this

#Adding adjust flag ##

system("./plink2 --pfile Fcy_rs_LD --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --adjust --out Fc_LD")

#Notes 26/01/22 ####
#So far these results are randomly pruned for linkage disequilibrium, could miss out on the more significant SNPs see snp_annotating script and haploview for extensions
#For now let's visualize our top snps with no pruning/prioritisation

#SNP Case-Control Plotting ####

#Plot top 'protective' SNP and top suscepitible SNP
#Case control analysis 

#Steps ####

#Update genotype information to RsID 
#see pos2rsID script
#Use Fcy_rs pgen files
#Get genotype case control info from --hardy
#Filter for sig SNPs
#Recode poss
#Make stacked bar chart
ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  geom_bar(stat="identity")

#Plink2 version hardy ####

system("./plink2 --pfile Fcy_rs --pheno pheno_update.txt --pheno-name Diagnosed --hardy --keep-if Diagnosed == 1 --out nTD")#139 controls

nTD <- fread("nTD.hardy")

#use keep if to make a diagnosed table, can then compare to undiagnosed table

system("./plink2 --pfile Fcy_rs --pheno pheno_update.txt --pheno-name Diagnosed --hardy --keep-if Diagnosed == 2 --out TD") #173 cases

TD <- fread("TD.hardy")

#27/01/22 Notes ####
# Updating script to stick with PLINK2 and use the Fcy_rs files #Also deleting binary code conversion matrix
# Assoc counts/frequency option don't work either, as the output does not split by genotype.
#The only option I can find that does is --hardy. SNPtest might also be better suited

#Merging TD nTD ####
#Need to reformat into one dataframe that we can plot 
#Merging will add .x .y suffix to colnames, add suffix option of .case for example when merging dataframes
#might need to split into separate observations so same sample repeated then a code column for case/control wide >long

#remove extra cols for now
nTD <- nTD[,2:7]
TD <- TD[,2:7]

#Merge hardy results ####
genotype <- left_join(nTD, TD, by = "ID", suffix = c(".cont",".case"))

#split into long format
#Actually because I hate reformatting adding column Pheno to both dataframes, then binding

pheno <- rep("Cont",487)
nTD$pheno <- pheno
pheno <- rep("Case",487)
TD$pheno <- pheno

genotype <- bind_rows(nTD, TD)
#double check it keeps same ref alt allele in each

#Plotting ####

#Current, most sig 'protective' SNP
#rs11590932 OR 0.4

#y axis we would want to be counts format to long again
geno <- melt(genotype, id.vars = c("ID", "pheno"), measure.vars = c(4:6))

colnames(geno) <- c("ID", "pheno", "genotype", "count")
geno$pheno <- as.factor(geno$pheno)

#Protective SNPs Plots ####
geno %>% filter(ID == "rs11590932") %>% ggplot(aes(x=genotype, y=counts, fill=pheno)) + geom_bar(stat = "identity") + theme_light()
#Having 1 copy of the allele het higher in controls than cases
#Same for two copies of the allele

geno %>% filter(ID == "rs180978155") %>% ggplot(aes(x=genotype, y=counts, fill=pheno)) + geom_bar(stat = "identity") + theme_light()


#Susceptibility SNP plots ####
geno %>% filter(ID == "rs12040409") %>% ggplot(aes(x=genotype, y=counts, fill=pheno)) + geom_bar(stat = "identity") + theme_light()
#Having Ax (minor) allele is much more associated with case of typhoid fever
#You'd want first column to be somewhat even
#then looking at the effects as the SNP occurs either one allele (het/two alleles) that increases the likelihood of case of typhoid fever

#Also test rs115297732, rs10917740, rs12040409, rs796681563
#these are in LD with eachother I think

#Polygenic Risk score, haplotypes, FcAR still to go...

#Just Typhoid ####

system("./plink2 --pfile Fcy_rs --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5, --remove-if Challenge == 3 --covar-variance-standardize --glm --out fcy_ty")

gwasResults <- fread("fcy_ty.Diagnosed.glm.logistic.hybrid")
#Filter TEST
gwasResults$TEST <- as.factor(gwasResults$TEST)
gwasResults <- filter(gwasResults, TEST == "ADD")
sig <- filter(gwasResults, P <0.05)

#rs115297732 #susceptible ##
geno %>% filter(ID == "rs115297732") %>% ggplot(aes(x=genotype, y=counts, fill=pheno)) + geom_bar(stat = "identity") + theme_light()

#rs11590932
geno %>% filter(ID == "rs11590932") %>% ggplot(aes(x=genotype, y=counts, fill=pheno)) + geom_bar(stat = "identity") + theme_light()

levels(geno$genotype) <- c("0", "1", "2")
levels(geno$pheno) <- c("TD", "nTD")

#dodged graph ##
geno %>% filter(ID == "rs115297732") %>% ggplot(aes(x=genotype, y=counts, fill=pheno)) + geom_bar(stat="identity", color="black", position=position_dodge()) + theme_light()

# eqtl

setwd("~/GWAS_22/Fc_receptor/data")


pheno_exprs <- fread("pheno_express_T1T2.txt")
geno <- fread("Fc_genoD0.txt")




