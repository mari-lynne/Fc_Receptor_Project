#! bin/bash'

#file is merged from T1P1 VAST_PATCH imputation runs
#cd Plink
# --double-id #make a dummy ID list later to figure out the dupicates

./plink2 --vcf all_enteric_bi.vcf.gz --const-fid --out enteric

#see R script post-impute for Plink commands
#update IIDs, added 26 to duplicates in P1/T1

#QC

./plink2 --pfile enteric --geno 0.01 --maf 0.05 --hwe 0.000001 --make-pgen --out enteric_QC

./plink2 --pfile enteric_QC --mind 0.01 --make-pgen --out enteric_QC

#./plink2 --pfile enteric_QC --pheno pheno_rn.txt --pheno-name Diagnosed --covar covar_rn.txt --covar-name Sex, Challenge, Dose, Study, Vaccine, Age --pca

./plink2 --pfile enteric_QC --pheno pheno_update.txt --pheno-name Diagnosed --covar covar_update.txt --covar-name Sex, Challenge, Dose, Vaccine, Age --pca

./plink2 --pfile enteric_QC --pheno pheno_update.txt --pheno-name Diagnosed --covar covar_pca2.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --glm  --covar-variance-standardize --adjust --out enteric_assoc

#replace categorical NAs with NONE
#make all categorical
#Replace numeric NAs with NA #might be getting involved
#replace FALSE in PC cols with NA #see r script for more troubleshooting, no joy try in plink 1.9

./plink2 --pfile enteric_QC --make-bed -out enteric_QC

./plink --bfile enteric_QC --pheno pheno2.txt --pheno-name Diagnosed --allow-no-sex --covar covar2.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --logistic --adjust --out enteric_assoc

