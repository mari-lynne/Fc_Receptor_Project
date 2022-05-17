#load packages ####
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(tidylog)
library(limma)
library(stringi)
library(janitor)
library(stringr)
library(splines)


#Just plot FCR splines models
setwd("~/RNA/Daniel")
#Data generated in time_course.R script
save.image(file= "Fcsplines.RData")
load("Fcsplines.RData")

#When generalising script, call function to count ncols in dfs - save as numeric var (floating point?)

#FCR data set up ####

#Get FCR data
pData2 <- pheno_exprs[,c(1:14)]
exprs2 <- pheno_exprs[,c(14:10461)]
 #FIX: removed time_ID col
#exprs2 <- exprs2 %>% mutate_at(c(2:10461), a)s.numeric

Fc_exprs <- exprs2 %>%
                select(fcar,fcgr1b, fcgr2a, fcgr2a_2, fcgr3a, fcrla, fcrlb)
  Fc_exprs <- sapply(Fc_exprs, as.numeric) %>%
                    cbind(pData2)

#Plot raw data across time ####
  #Focus on first two weeks after challenge
table(Fc_exprs$time)
Fc_exprs %>% ggplot(aes(x=time, y=fcar, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 15) +
  scale_colour_manual(values=c("seagreen3", "sandybrown"))


#Other plots to make = volcano plots and spline regression 

#Spline regression ####

#Spline Plots ####
Fc_exprs %>% ggplot(aes(x=time, y=fcgr1b, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 5)) +xlab("Time (days)") + ylab(expression(paste("FC",gamma,"R1",beta, " Expression")))

Fc_exprs %>% ggplot(aes(x=time, y=fcgr2a, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14) + ylim(10,13.3) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 3)) +
  xlab("Time (days)") + ylab(expression(paste("Fc",gamma,"R2",alpha, " Expression")))

Fc_exprs %>% ggplot(aes(x=time, y=fcgr3a, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 4)) + 
  xlab("Time (days)") + ylab(expression(paste("Fc",gamma,"R3",alpha, " Expression"))
#Not much interesting with 3a


Fc_exprs %>% ggplot(aes(x=time, y=fcar, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 4)) +
  xlab("Time (days)") + ylab(expression(paste("Fc",alpha,"R", " Expression")))


Fc_exprs %>% ggplot(aes(x=time, y=fcrla, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 3)) + 
  xlab("Time (days)") +
  ylab(expression(paste("FcRL",alpha, " Expression")))


Fc_exprs %>% ggplot(aes(x=time, y=fcrlb, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 3)) + 
  xlab("Time (days)") + ylab(expression(paste("FcRL",beta, " Expression")))








#Steps
#Make design matrix
#Set up splines model
#Fit model
#Evaluate model
#Plot splines regression

#We want to run two separate splines models, one for TD partipants across time and another for nTD
#Remove later time points for now

#Reform eset obeject so we can input to limma easily:
#The number of rows in phenoData must match the number of columns in assayData
#Row names of phenoData must match column names of the matrix / matrices in assayData.

#Resplt Fc_exprs
pData2 <- Fc_exprs[,c(8:21)]
exprs2 <- Fc_exprs[,c(1:7)]

#Make Time_ID row names in pheno data
rownames(pData2) <- pData2$time_id 

#Re Transpose exprs data frame 
exprs2 <- transpose(exprs2)
names(exprs2) <- exprs2[1,] %>% as.vector(mode = "character")
exprs2 <- exprs2[-c(1),] 
#recode to numeric
exprs2 <- sapply(exprs2, as.numeric)


#Make design matrix ####
#can be done manually by selecting appropriate vars as a df
#or using model.matrix function
#Check matrix against design on paper 

#Set up vars for design matrix ####
#make time as numeric
pData2$time <- as.numeric(pData2$time)
pData2$study_arm <- as.factor(pData2$study_arm)
#make new vaccine factor (will give more df)
pData2$vax <- ifelse(pData2$study_arm == "T1",0,
                     ifelse(pData2$study_arm == "Placebo",0,
                            ifelse(pData2$study_arm == "M01ZH09", 1,
                                   ifelse(pData2$study_arm == "Ty21a", 1,NA))))
pData2$vax <- as.factor(pData2$vax)

design_full <- model.matrix(~diagnosis + time + age + sex + vax, data = pData2)

#Splines ####
library(splines)


#We also need to add splines points to our design matrix so lm function can use these to estimate the model
#Make a separate splines matirx or df before combining this with our original design

design_splines <- ns(pData2$time, df = 5, intercept = TRUE)
#Boundary.knots = range(x))


#Full matrix ####
#design_splines_full <- model.matrix(~design_full + design_splines)
#Just bind matricies as they're compatible
design_splines_full <- cbind(design_full, design_splines)

#Model per diagnosis across time using splines ####
library(limma)
lm <- lmFit(exprs2, design = design_splines)

summary(lm)

