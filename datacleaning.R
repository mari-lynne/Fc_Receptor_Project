setwd("C:/Users/mjohnson/Documents/genetics/plink-1.07-dos")

library("dplyr")

#load data
p1master <- read.csv("p1tyger_pheno.csv", fileEncoding = "UTF-8-BOM")
p1extra <- read.csv("p1_clean.csv" , fileEncoding = "UTF-8-BOM")


#recode sample IDs

p1master$Outcome <- as.factor(p1master$Outcome)
#1master$Sample_ID <- as.factor(p1master$Sample_ID)#
str(p1master)


p1extra$Outcome <- as.factor(p1extra$Outcome)
#1extra$Sample_ID <- as.factor(p1extra$Sample_ID)
str(p1extra)


#recode typhoid diagnosis

p1master$Outcome <- ifelse(p1master[,2] == "Typhoid", 1, ifelse(p1master[,2] == "No Typhoid", 0,NA))

#merge table (not very good)

newtable <- left_join(p1master, p1extra, by = "Sample_ID")

write.csv(newtable, "newtable.csv")

#change sexinto 

newtable2 <- read.csv("newtable2.csv",  fileEncoding = "UTF-8-BOM")


newtable2$Sex <- ifelse(newtable2$Sex == "M", 1, ifelse(newtable2$Sex == "F", 2, NA))

str(newtable2)

write.csv(newtable2, "newtable2.csv")

#merge so just with genotype study samples

gensamples <- read.csv("samplesp1tyger.csv",  fileEncoding = "UTF-8-BOM")

str(gensamples)


gentable <- left_join(gensamples, newtable2, by = "Sample_ID")

write.csv(gentable, "gentable.csv")


gentable2 <- read.csv("gentable2.csv",  fileEncoding = "UTF-8-BOM")

gentable2$Outcome <- ifelse(gentable$Outcome == "1", 2, ifelse(gentable$Outcome == "0", 1, NA))

write.csv(gentable2, "gentable2.csv")


#load t1t2 data

t1t2key <- read.csv("sample_key_t1t2.csv", fileEncoding = "UTF-8-BOM")

t1t2pheno <- read.csv("t1t2phenotype.csv", fileEncoding = "UTF-8-BOM")

t1t2new <- left_join(t1t2key, t1t2pheno, by = "LabID")

#recode variables

t1t2new$Outcome <- ifelse(t1t2new$Outcome == "No Typhoid", 1, ifelse(t1t2new$Outcome == "Typhoid", 2, NA))

t1t2new$Sex <- ifelse(t1t2new$Sex == "M", 1, ifelse(t1t2new$Sex  == "F", 2, NA))

View(t1t2new)

write.csv(t1t2new, "t1t2new.csv")

#more recoding

t1cov <- read.csv("t1t2new.csv" , fileEncoding = "UTF-8-BOM")

t1cov$Vaccine <- ifelse(t1cov$Vaccine == "Placebo", 0, ifelse(t1cov$Vaccine == "M01ZH09", 1, ifelse(t1cov$Vaccine == "Ty21a", 2, NA)))


t1cov$Strain.or.dose <- ifelse(t1cov$Strain.or.dose == "Quailes 1-5x10^4 ", 4, ifelse(t1cov$Strain.or.dose == "Quailes 1-5x10^3", 5, NA))

#study t1 = 3, t2 =4

write.csv(t1cov, "t1cov2.csv")


#recode covariate dose groups p1 tyger

p1cov <- read.csv("p1tygercovariate.csv" , fileEncoding = "UTF-8-BOM")

p1cov$Strain.or.dose <- ifelse(p1cov$Strain.or.dose == "Paratyphi A 1-5x10^3 ", 2, ifelse(p1cov$Strain.or.dose == "Paratyphi A 0.5-1x10^3", 1, ifelse(p1cov$Strain.or.dose == "TN", 3, ifelse(p1cov$Strain.or.dose == "Quailes 1-5x10^4 ", 4, NA))))


p1cov$Study <- ifelse(p1cov$Study == "P1", 1,ifelse(p1cov$Study == "Tyger", 2, NA))


write.csv(p1cov, "p1cov.csv")




rename.values(p1cov, Paratyphi A 1-5x10^3="Phigh")

install.packages("plyr")
library(plyr)
str(p1cov)

recode(p1cov$Strain.or.dose, Paratyphi A 1-5x10^3 ="Phigh")

mapvalues(x, from = c("beta", "gamma"), to = c("two", "three"))

p1cov$Strain.or.dose <- mapvalues(p1cov$Strain.or.dose, from = c("Paratyphi A 1-5x10^3", "Paratyphi A 0.5-1x10^3"), to = c("Phigh", "Plow"))








