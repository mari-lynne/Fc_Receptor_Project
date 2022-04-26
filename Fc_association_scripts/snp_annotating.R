#dSNP ####

#identify candidate genes by annotation, then LD

dsnp <- fread("rs_list_dsnp")

dsnp$hg38.snp151.class <- as.factor(dsnp$hg38.snp151.class)
dsnp$hg38.snp151.func <- as.factor(dsnp$hg38.snp151.func)

levels(dsnp$hg38.snp151.func) 

dsnp <- dsnp %>% mutate(dup_check = stri_duplicated(dsnp$hg38.snp151.chromEnd))
dsnp <- dsnp %>% filter(dup_check == "FALSE")

dsnp <- rename(dsnp,SNP = hg38.snp151.name)

annotate <- left_join(sig, dsnp, by = "SNP")

#get lead snps (basically interesting ones)
annotate$hg38.snp151.func <- as.factor(annotate$hg38.snp151.func)

lead_snps <- annotate %>% filter(hg38.snp151.func != "unknown", hg38.snp151.func != "intron")

lead_snps <- annotate %>% filter(str_detect(hg38.snp151.func, "near|coding"))

#just get rsID :rsID of the lead SNPs
#chr : chromosome
#pos : genomic position (hg19) #liftover

lead_snps <- select(lead_snps, SNP, CHR, BP)

write.table(lead_snps, file = "lead_snps.txt", sep = "\t", quote = F, col.names = F, row.names = F)

#Fuma updates 

#The pipeline currently supports human genome hg19. If your input file is not based on hg19, please update the genomic position using liftOver from UCSC. However, there is an option for you!! When you provide only rsID without chromosome and genomic position, FUMA will extract them from dbSNP build 146 based on hg19. To do this, remove columns of chromosome and genomic position or rename headers to ignore those columns. Note that extracting chromosome and genomic position will take extra time.

fuma_snps <- select(fuma, -CHR, -BP)
write.table(fuma_snps, file = "fuma_snps.txt", sep = "\t", quote = F, col.names = F, row.names = F)

#Liftover lead
# BED= chr:Start:End chr4 100000 100001

#LD.annotate ####

#where "PathToSnpFiles" is the path to the folder containing all data file
#"annot.gff3" is the file containing genomic coordinates and annotations for features (most often genes)
#"candidate" is a list of chromosomes and positions for candidate polymorphisms
#"type" is the feature (mRNA, CDS, gene)
#"thr" is the threshold for r2, "output" is an output name specified by the user, an
#"SNP_Map" is a map file indicating chromosome and positions for all SNPs genotyped using the SNP-chip.
#genomic coordinates 


#Mari Version #####

#First independent filter R2 = <0.6
system("./plink2 --pfile Fcy_rs --indep-pairwise 50 5 0.6 --out Fcy_ind1")

#374/487 variants removed = 113 use these to annotate
system("./plink2 --pfile Fcy_rs --extract Fcy_ind1.prune.in --make-pgen --out Fcy_ind1")

#make annotation file for https://www.snp-nexus.org/v4/guide/

pvar <- fread("Fcy_ind1.pvar")
#dsnp rsID

nexus_snp <- nexus[,3]

nexus_snp <- mutate(nexus_snp, dbsnp = rep("dbsnp"))
nexus_snp <- select(nexus_snp, dbsnp, ID)

write.table(nexus_snp, file = "nexus_snp.txt", sep = "\t", quote = F, col.names = F, row.names = F)

#output ####

#could filter from CAD?
#LD musings ####
#we have the filter snps that are in <0.6 in hapmap/1000genome #then proceed
#straight up filter snps <0.2 (maybe increase to 0.3)
#PHRED of <12 = likely to be deleterious
#(?<=chr1:)\\d+

coding_exon <- coding %>% filter(str_detect(Predicted_Function, "^((non-coding).)*"))

coding_exon <- coding %>% filter(str_detect(Predicted_Function,"^((?!non-coding).)*$"))

coding_exon <- coding_exon %>% filter(str_detect(Predicted_Function,"^((?!intronic).)*$"))

#137 SNPs coding filtered R2 of <0.6
write.table(coding_exon$ID, file = "coding_exon.txt", sep = "\t", quote = F, col.names = F, row.names = F)

#use this in extract to perform assoc with Fcy_Rs #test_GWAS in assoc ##

#No LD, Just filter on function ####
#Then check for LD ##

#upload the 400 snps to nexus then redo coding/promoter filter

system("./plink2 --pfile Fcy_rs --indep-pairwise 50 5 0.8 --out Fcy_ind2")

#374/487 variants removed = 113 use these to annotate #0.8 gives 143 varients
#I think include them all to have our best chance of coding, then filter 0.6
#make annotation file for https://www.snp-nexus.org/v4/guide/

#No LD
pvar <- fread("Fcy_rs.pvar")
#dsnp rsID

nexus_snp <- pvar[,3]

nexus_snp <- mutate(nexus_snp, dbsnp = rep("dbsnp"))
nexus_snp <- select(nexus_snp, dbsnp, ID)

write.table(nexus_snp, file = "nexus_snp_noLD.txt", sep = "\t", quote = F, col.names = F, row.names = F)

#maybe also keep snps in regulatory elements 

#rs1801274 SNP of interest, susceptibility to infection

#poss look into this 
https://github.com/molgenis/systemsgenetics/wiki/Downstreamer

