#PCA ####

#Recode Variables ####

#Numeric recoding
#Must use NA

#NONE = 0
#Ty21a = 1
#M01ZH09= 2
#Vi-PS = 3
#Vi-TT= 4

covar <- fread("covar_update.txt")

covar$Vaccine <- ifelse(covar$Vaccine == "None",0,
                            ifelse(covar$Vaccine == "Ty21a",2,
                                   ifelse(covar$Vaccine == "M01ZH09", 2,
                                          ifelse(covar$Vaccine == "Vi-PS", 3,
                                                 ifelse(covar$Vaccine == "Vi-TT", 4,NA)))))

#Recode strain

covar$Challenge <- as.factor(covar$Challenge)

levels(covar$Challenge)

#C, P = 3, T = 1, TN = 2
levels(covar$Challenge) <- c("0", "3", "1", "2")


#sex M=1, F=2
covar$Sex <- ifelse(covar$Sex == "M", 1, ifelse(covar$Sex == "F", 2, NA))

#Diagnosed 0>1 = undiagnosed, 1>2 = diagnosed
pheno$Diagnosed <- ifelse(pheno$Diagnosed == "0", 1, ifelse(pheno$Diagnosed == "1", 2, NA))

#update tables
write.table(covar, file = "covar_update.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote=F)

write.table(pheno, file = "pheno_update.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote=F)


#QC data in plink ####
./plink2 --pfile enteric --geno 0.01 --maf 0.05 --hwe 0.000001 --make-pgen --out enteric_QC

./plink2 --pfile enteric_QC --mind 0.01 --make-pgen --out enteric_QC

20476788 variants loaded from enteric.pvar.
Note: No phenotype data present.
Calculating allele frequencies... done.
--geno: 9228600 variants removed due to missing genotype data.
--hwe: 135 variants removed due to Hardy-Weinberg exact test (founders only).
5391382 variants removed due to allele frequency threshold(s)
(--maf/--max-maf/--mac/--max-mac).

#Number of SNPs remaining = 585,6671
#9228600 removed due to missing genotype data = imputation build/annotation issue #I think check MAF could be just rare imputed alleles

#crack on for now :/ see notebook for plans
#PCA ####
system("./plink2 --pfile enteric_QC --pheno pheno_update.txt --pheno-name Diagnosed --covar covar_update.txt --covar-name Sex, Challenge, Dose, Vaccine, Age --pca")

pca <- fread("./plink2.eigenvec")
eigenval <- scan("./plink2.eigenval")
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

#colour code pca by sex, ethnicity
rm(covar_pca)
pca_covar <- left_join(pca, covar, by = "IID") #has the extra RPT samples
missing_pca <- anti_join(covar, pca, by = "IID") #791

write.table(pca_covar, file = "pca_covar.txt", row.names = FALSE, col.names = TRUE, sep="\t", quote=F)


