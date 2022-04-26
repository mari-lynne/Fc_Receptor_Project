##Mari expression attemts#####

#make first column a column
logFC_12 <- setDT(logFC_12, keep.rownames = TRUE)[]
setDT(logFC_24, keep.rownames = TRUE)[]

library(dplyr)


#Find Fc expression data 12h 

FCGR_expression <-filter(logFC_12, grepl('FCG', rn))
FCRL_expression <-filter(logFC_12, grepl('FCR', rn))
FCAR_expression <-filter(logFC_12, grepl('FC', rn))

test <- ggplot(Fcgr_expression,aes(rn))

test + geom_bar() 

write.csv(Fcgr_expression, "FCGR_12hFC.csv")

getwd()

setwd("C:/Users/mjohnson/Documents/RNA")

FCG_12FC <-read.csv("FCGR_12hFC.csv")

test <- ggplot(FCG_12FC,aes(Participant_ID))

FCG_12FC %>% filter("FCGR2A")

test + geom_bar() 
str(FCG_12FC)

library(ggplot2)

FCG_12FC %>% ggplot(aes(x=Participant_ID, y= FCGR2A, color = Participant_ID)) +geom_bar(stat= "identity") +theme(legend.position = "none") +xlab("Participants") +ylab("LogFC expression")+ ggtitle("FCGR2A")

FCGR2A


FCGR2B <- FCG_12FC %>% ggplot(aes(x=Participant_ID, y= FCGR2B, color = Participant_ID)) +geom_bar(stat= "identity") +theme(legend.position = "none") +xlab("Participants") +ylab("LogFC expression")+ ggtitle("FCGR2B")

FCGR2B


FCGR3A<- FCG_12FC %>% ggplot(aes(x=Participant_ID, y= FCGR2B, color = Participant_ID)) +geom_bar(stat= "identity") +theme(legend.position = "none") +xlab("Participants") +ylab("LogFC expression")+ ggtitle("FCGR3A")

FCGR3A

FCGR3A<- FCG_12FC %>% ggplot(aes(x=Participant_ID, y= FCGR2B, color = Participant_ID)) +geom_bar(stat= "identity") +theme(legend.position = "none") +xlab("Participants") +ylab("LogFC expression")+ ggtitle("FCGR3A")

FCGR3A

library(stringr)

#match with participant data and outcomes

FCG_12FC$Participant_ID <- FCG_12FC$Participant_ID %>% str_remove_all("[P]")


meta <- read.csv("t1t2new.csv")

str(meta)

meta$Participant_ID <- as.character(meta$Participant_ID)

FCG_12FC_annotated <- FCG_12FC %>% inner_join(meta, by = c("Participant_ID"))

FCG_12FC_annotated$Outcome <- as.factor(FCG_12FC_annotated$Outcome)

FCGR3A_outcome <- FCG_12FC_annotated %>% ggplot(aes(x=Participant_ID, y= FCGR3A, color = Outcome)) +geom_bar(stat= "identity") +xlab("Participants") +ylab("LogFC expression")+ ggtitle("FCGR3A")

FCGR3A_outcome

FCG_12FC_annotated %>% ggplot(aes(x=Participant_ID, y= FCGR2A, color = Outcome)) +geom_bar(stat= "identity") +xlab("Participants") +ylab("LogFC expression")+ ggtitle("FCGR2A")



FCGR2B_outcome <- FCG_12FC_annotated %>% ggplot(aes(x=Participant_ID, y= FCGR2B, color = Outcome)) +geom_bar(stat= "identity") +xlab("Participants") +ylab("LogFC expression")+ ggtitle("FCGR2B")

FCGR2B_outcome


#cibersort

install.packages("remotes")
remotes::install_github("icbi-lab/immunedeconv")

library(immunedeconv)

#load gene expression matrix in

gene_exprs <- read.table("CiberSortInput.txt", header = TRUE, row.names = TRUE)



out <- immunedeconv::deconvolute(gene_exprs, "epic")



install.packages("ADAPTS")
library(ADAPTS)

data("LM22")
