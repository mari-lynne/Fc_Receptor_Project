#FcaR testing

#Get FcaR coordinates GChr38
#filter enteric_QC


#54874231..54891420 #minus/plus 10kb

system("./plink2 --pfile enteric_QC --chr 19 --from-bp 54864231 --to-bp 54901420 --make-pgen --out fcar") #191 varients

system("./plink2 --pfile fcar --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4 --covar-variance-standardize --glm --out fcar")

gwasResults <- fread("fcar.Diagnosed.glm.logistic.hybrid")
#Filter TEST
gwasResults$TEST <- as.factor(gwasResults$TEST)
gwasResults <- filter(gwasResults, TEST == "ADD")
sig <- filter(gwasResults, P <0.05)

#one snp is p = 0.06