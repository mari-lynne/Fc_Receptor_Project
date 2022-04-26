library(qqman)

#Logistic regression ###

system("./plink2 --pfile enteric_QC --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --out enteric_assoc")

#Results written to enteric_assoc.Diagnosed.glm.logistic.hybrid

gwasResults <- fread("enteric_assoc.Diagnosed.glm.logistic.hybrid")

#Filter TEST
gwasResults$TEST <- as.factor(gwasResults$TEST)

gwasResults <- filter(gwasResults, TEST == "ADD")

#Plot
colnames(gwasResults) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1","FIRTH?","TEST","OBS_CT","OR", "LOG(OR)_SE","Z_STAT","P","ERRCODE")

manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 8), cex = 0.6, cex.axis = 0.9, col = c("#34b368", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c(1:20, "P", "Q"))

#sort by P value
gwasResults
view(gwasResults)

#SNP annotation still dodgy, is this from indexing?
#VCF files convert using dSNP from usc track

#Fcy Receptor associations #CHR38!

system("./plink2 --pfile enteric_QC --chr 1 --from-bp 161505430 --to-bp 161678022 --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --out enteric_assoc")

gwasResults <- fread("enteric_assoc.Diagnosed.glm.logistic.hybrid")

#Filter TEST
gwasResults$TEST <- as.factor(gwasResults$TEST)

gwasResults <- filter(gwasResults, TEST == "ADD")

#Plot
colnames(gwasResults) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1","FIRTH?","TEST","OBS_CT","OR", "LOG(OR)_SE","Z_STAT","P","ERRCODE")

#sort by P value
sig <- filter(gwasResults, P <0.05)

#filter LD ####

system("./plink2 --pfile enteric_QC --chr 1 --from-bp 161505430 --to-bp 161678022 --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --out enteric_assoc")

system("./plink2 --pfile enteric_QC --indep-pairwise 50 5 0.2 --out enteric_LD")

system("./plink2 --pfile enteric_QC --extract enteric_LD.prune.in --make-pgen --out enteric_LD")

#FC assoc with LD ####
#currently 172kb
#add 100kb either side 
#100kb less, then plus

#100kb = 100,000

system("./plink2 --pfile enteric_LD --chr 1 --from-bp 160505430 --to-bp 162678022 --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --out Fc_LD")

#That gives 261 varients
#1KB extra either side

#Fcy updated co-ordinates glm ####

#from viewer 
#minus/plus 10kb

#FCGR2a = 161504430
#FCRLB = 161729143
#=currently a 224kb window


system("./plink2 --pfile enteric_LD --chr 1 --from-bp 161504430 --to-bp 161729143 --make-pgen --out Fcy")

system("./plink2 --pfile Fcy --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --out Fcy_LD")

#26 variants LD

gwasResults <- fread("Fcy_LD.Diagnosed.glm.logistic.hybrid")

#Filter TEST
gwasResults$TEST <- as.factor(gwasResults$TEST)

gwasResults <- filter(gwasResults, TEST == "ADD")
sig <- filter(gwasResults, P <0.05)

#4 significant snps
#161704842 near FCRLA

#No LD pruning Fcy - Fcy_2 ####
#see FUMA tool

system("./plink2 --pfile enteric_QC --chr 1 --from-bp 161504430 --to-bp 161729143 --make-pgen --out Fcy_2")

#487 varients
system( )

gwasResults <- fread("Fcy.Diagnosed.glm.logistic.hybrid")

#Filter TEST
gwasResults$TEST <- as.factor(gwasResults$TEST)

gwasResults <- filter(gwasResults, TEST == "ADD")
sig <- filter(gwasResults, P <0.05)

#72 significant SNPs
#see pos2rsID for updating SNP IDs

#Fcy_rs is the updated file :)

#redo assoc for FUMA input ####

system("./plink2 --pfile Fcy_rs --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --out Fcy_rs")

#check glm name 

gwasResults <- fread("Fcy_rs.Diagnosed.glm.logistic.hybrid")

#Filter TEST
gwasResults$TEST <- as.factor(gwasResults$TEST)

gwasResults <- filter(gwasResults, TEST == "ADD")

sig <- filter(gwasResults, P <0.05)


#sort col names for FUMA
#SNP | snpid | markername | rsID: rsID
#CHR | chromosome | chrom: chromosome
#BP | pos | position: genomic position (hg19)
#A1 | effect_allele | allele1 | alleleB: affected allele (ALT)
#A2 | non_effect_allele | allele2 | alleleA: another allele (REF)
#P | pvalue | p-value | p_value | pval: P-value (Mandatory)
#OR: Odds Ratio
#Beta | be: Beta #log_SE
#SE: Sta

colnames(sig) <- c("CHR", "BP", "SNP", "REF", "A2", "A1","FIRTH?","TEST","OBS_CT","OR", "Beta","Z_STAT","P","ERRCODE")

fuma <- select(gwasResults, SNP, CHR, BP, A1, A2, P, OR, Beta)

write.table(fuma, file = "fuma.txt", sep = "\t", quote = F, col.names = T, row.names = F)

fuma_sig <- select(sig, SNP, CHR, BP, A1, A2, P, OR, Beta)

write.table(fuma_sig, file = "fuma_sig.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#If fuma doesn't work download genseq to include SNPs in interesting regions 

#Need to predefine lead SNPs
#get genome track

#dSNP ####

#identify candidate genes by annotation, then LD
#see snp_annotating script
system("./plink2 --pfile Fcy_rs --extract coding_exon.txt --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --covar-variance-standardize --glm --out Fcy_interest")

#30 coding Varients

gwasResults <- fread("Fcy_interest.Diagnosed.glm.logistic.hybrid")

#Filter TEST
gwasResults$TEST <- as.factor(gwasResults$TEST)
gwasResults <- filter(gwasResults, TEST == "ADD")
sig <- filter(gwasResults, P <0.05)

#leaves 4 rows remaining 
#test these for LD
#In a candidate gene study, trying to find snps that have an effect not just a casual variant

#extract IDs, compare to coding annotate

sig <- select(sig, ID, REF, ALT, OR, P)
sigsnps <- semi_join(sig, coding_exon, by = "ID")


#Test without Vi-PS, Vi-TT

--remove-if <phenotype/covariate name> <operator> <value>
--remove-if Vaccine = 4
system("./plink2 --pfile Fcy_rs --extract coding_exon.txt --pheno pheno_update.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --remove-if Vaccine = 4 --glm --out Fcy_noTT")

#10 covars in glm more than 10th of control

#make new pheno file
pheno <- fread("pheno_update.txt")
covar <- fread("pca_covar.txt")

covar <- filter(covar, Vaccine != "3" & Vaccine != "4")
#removed 80 samples, sounds about right
#239 rows left
#I also have added PATCH data left to work out

pheno <- left_join(covar, pheno, by = "IID")

pheno <- select(pheno, IID, Diagnosed)

write.table(pheno, file = "noVi_pheno.txt", sep = "\t", quote = F, col.names = T, row.names = F)
#rerun analysis ####  

system("./plink2 --pfile Fcy_rs --extract coding_exon.txt --pheno noVi_pheno.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --remove-if Vaccine = 4 --glm --covar-variance-standardize --out Fcy_noTT")

#No LD, Just filter on function ####
#Then check for LD ##


system("./plink2 --pfile Fcy_rs --extract coding_exon.txt --pheno noVi_pheno.txt --pheno-name Diagnosed --covar pca_covar.txt --covar-name Sex, Challenge, Dose, Vaccine, Age, PC1, PC2, PC3, PC4, PC5 --remove-if Vaccine = 4 --glm --covar-variance-standardize --out Fcy_noTT")

