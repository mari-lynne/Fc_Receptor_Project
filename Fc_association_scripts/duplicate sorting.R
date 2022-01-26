#Check for duplicates and update

#Load files
#Update pheno ####
pheno <- fread("pheno.txt")
covar <- fread("covar.txt")
psam <- fread("enteric.psam") 
#psam <- fread("enteric_QC.psam") dont need to do

#check for duplicates 

U_P = unique(pheno$GENOTYPING_ID) #309
U_C = unique(covar$GENOTYPING_ID) #309

pheno_covar <- bind_cols(pheno, covar) #match up by rows

pheno_covar <- select(pheno_covar, -GENOTYPING_ID...5)

#basically don't change the genotyping ID u idiot
    
#reorder by factor to match psam ####
#then update duplicates
pheno_covar$Study <- as.factor(pheno_covar$Study)

levels(pheno_covar$Study) #not order

category_order <-c("P1", "TYGER", "T1", "T2", "VAST", "PATCH")

pheno_covar <- pheno_covar %>% 
  arrange(factor(Study, levels = category_order))

#Update duplicate IDs ####

#edit IID column
pheno_covar <- pheno_covar %>% rename(IID = GENOTYPING_ID...2)

pheno_covar <- pheno_covar %>% mutate(dup_check = stri_duplicated(pheno_covar$IID))
#highlights first duplicate entry

pheno_covar$IID <- ifelse(pheno_covar$dup_check == TRUE, 
                     paste(pheno_covar$IID , "26",sep = ""),
                     paste(pheno_covar$IID, "",sep = ""))

upc <-unique(pheno_covar$IID) #318 unique values :)

#Make pheno/covar tables
pheno <- select(pheno_covar, IID, Diagnosed)
covar <- select(pheno_covar, -Diagnosed)

#update psam file ####
#(remove hashtag in txt)
#remove first part of string before _

psam psam$IID <- sub(".+?_", "", psam$IID)

#add 26 to IID
psam <- psam %>% mutate(dup_check = stri_duplicated(psam$IID))
dups <- psam %>% filter(dup_check == "TRUE")

psam$IID <- ifelse(psam$dup_check == TRUE, 
                          paste(psam$IID , "26",sep = ""),
                          paste(psam$IID, "",sep = ""))

U_P <- unique(psam$IID) #319 woo

#check with left join

test <- left_join(psam, pheno, by = "IID")
missing1 <- anti_join(psam, pheno, by = "IID") #IDs not in pheno = RPT
missing2 <- anti_join(pheno, psam, by = "IID") #IDs not in psam 791
missing3 <- anti_join(psam, covar, by = "IID") #RPT
#Pretty much sorted :)

#Write files ####
psam <- select(psam, -dup_check)
covar <- select(covar, -dup_check)
pheno <- select(pheno, -dup_check)

write.table(psam, file = "enteric.psam", row.names = F, quote = F, sep = "\t")
write.table(covar, file = "covar_update.txt", row.names = F, quote = F, sep = "\t")
write.table(pheno, file = "pheno_update.txt", row.names = F, quote = F, sep = "\t")


