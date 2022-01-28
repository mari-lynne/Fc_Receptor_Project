#FcaR testing

#Get FcaR coordinates GChr38
#filter enteric_QC


#54874231..54891420 #minus/plus 10kb

system("./plink2 --pfile enteric_QC --chr 19 --from-bp 54864231 --to-bp 54901420 --make-pgen --out fcar") #191 varients

system("./plink2 --pfile fcar --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5, --covar-variance-standardize --glm --out fcar")

#just typhoid ####
system("./plink2 --pfile fcar --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5, --remove-if Challenge == 3 --covar-variance-standardize --glm --out fcar")

gwasResults <- fread("fcar.Diagnosed.glm.logistic.hybrid")
#Filter TEST
gwasResults$TEST <- as.factor(gwasResults$TEST)
gwasResults <- filter(gwasResults, TEST == "ADD")
sig <- filter(gwasResults, P <0.05)
#one snp is p = 0.06
#with just typhoid, p = 0.045

#SNP rs4560030 basically not present in afr or asian, but is in EUR 10% #Could that be indicative of selection pressure?

#Just Para ####
system("./plink2 --pfile fcar --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Dose, Age, PC1, PC2, PC3, PC4, PC5, --keep-if Challenge == 3 --covar-variance-standardize --glm --out fcar")
#try --remove-if Challenge == 3
#just typ
gwasResults <- fread("fcar.Diagnosed.glm.logistic.hybrid")
#Filter TEST
gwasResults$TEST <- as.factor(gwasResults$TEST)
gwasResults <- filter(gwasResults, TEST == "ADD")
sig <- filter(gwasResults, P <0.05)
#29 cases, 29 controls - need to get the rest of the paradata! they would be rechallenged in PATCH
