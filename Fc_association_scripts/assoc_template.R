#general assoc testing

system("./plink2 --pfile Fcy_rs --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --out Fcy_rs")

#30 Varients

gwasResults <- fread("Fcy_rs.Diagnosed.glm.logistic.hybrid")

#Filter TEST ###
gwasResults$TEST <- as.factor(gwasResults$TEST)
gwasResults <- filter(gwasResults, TEST == "ADD")
sig <- filter(gwasResults, P <0.05)

#just the one remaining snp with no vaccine :(
#0.500073, 0.004881 FCRLB

gwas.txt <- gwas$
  
