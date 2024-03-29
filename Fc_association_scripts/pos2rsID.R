#Update Chr position to rsID

#following advice from https://www.biostars.org/p/171557/

#Update rs_IDs ####

#which ones are interesting - use FUMA tool
#update rsIDs first

#see rs_list file use to update map
#chr:pos in column 1 and rs number in column 2.

rs<- fread("rs_list")
rs <- rs[,2:3]
write.table(rs, file="rs_list", sep = "\t", quote = F, col.names = F, row.names = F)

#update-map
#might have to update chr ID in pvar
pvar <- fread("Fcy_snps.pvar")

#edit ID col, split by :
test <- pvar[1:8,]
test$ID <- str_extract(test$ID, "(?<=chr1:)\\d+")

#could have also just copied the position column over daym
pvar$ID <- str_extract(pvar$ID, "(?<=chr1:)\\d+")

write.table(pvar, file = "Fcy_2.pvar", sep = "\t", quote = F, col.names = F, row.names = F) #add header back in manually

system("./plink2 --pfile Fcy_2 --update-name rs_list --make-pgen --out Fcy_rs")

#Remove duplicates ####
#extract unique snps from rs_list
rs <- rs %>% mutate(dup_check = stri_duplicated(rs$V1))
u_rs <- rs %>% filter(dup_check == "FALSE")

u_rs <- u_rs[,1:2]#unique snp names

write.table(u_rs, file = "u_rs.txt", sep = "\t", quote = F, col.names = F, row.names = F)

system("./plink2 --pfile Fcy_2 --update-name u_rs.txt --make-pgen --out Fcy_rs")

