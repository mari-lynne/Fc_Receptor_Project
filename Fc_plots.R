setwd(R)
#load packages ####
library(dplyr)
library(data.table)
library(tidyr)
library(tidylog)
library(limma)
library(stringi)
library(janitor)
library(stringr)
library(splines)
library(patchwork)

#Plots ####

#1) Volcano Plots
#2) FcR expression/splines plots
#3) FcR box plots 
#4) Genotype Plots
#5) Cis-EQTL plots

#2) FcR Splines ####

#Spline Plots ####
FC1B <- Fc_exprs %>% ggplot(aes(x=time, y=fcgr1b, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14.3) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 5)) +xlab("") + ggtitle(expression(paste("FC",gamma,"R1",beta))) + ylab("")

FC2A<- Fc_exprs %>% ggplot(aes(x=time, y=fcgr2a, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14.3) + ylim(10,13.3) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 3)) +
  xlab("") + ggtitle(expression(paste("FC",gamma,"R2",alpha))) + ylab("")

FC3A <- Fc_exprs %>% ggplot(aes(x=time, y=fcgr3a, colour = diagnosis)) +
  geom_point() + 
  xlim(0, 14.3) +
  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
  geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 4)) + 
  xlab("Time (days)") + ggtitle(expression(paste("FC",gamma,"R3",alpha))) + ylab("Expression\n")
#Not much interesting with 3a


                             
FCAR <- Fc_exprs %>% ggplot(aes(x=time, y=fcar, colour = diagnosis)) +
                               geom_point() + 
                               xlim(0, 14.3) +
                               scale_colour_manual(values=c("seagreen3", "sandybrown")) +
                               geom_smooth( method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 4)) + xlab("") + ggtitle(expression(paste("FC",alpha,"R"))) + ylab("Expression\n")


                             
FCRLA <- Fc_exprs %>% ggplot(aes(x=time, y=fcrla, colour = diagnosis)) +
                               geom_point() + 
                               xlim(0, 14.3) +
                               scale_colour_manual(values=c("seagreen3", "sandybrown")) +
                               geom_smooth( method = "lm",aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 3)) + 
                               xlab("Time (days)") + ggtitle(expression(paste("FCRL",alpha))) + ylab("")
                             
 FCRLB <- Fc_exprs %>% ggplot(aes(x=time, y=fcrlb, colour = diagnosis)) +
                  geom_point() + 
                  xlim(0, 14.3) +
                  scale_colour_manual(values=c("seagreen3", "sandybrown")) +
                  geom_smooth(method = "lm",  aes(group = diagnosis, colour = diagnosis), formula = y ~ splines::ns(x, df = 3)) + xlab("Time (days)") + ggtitle(expression(paste("FCRL",beta))) + ylab("")
  



List <- list(FCAR, FC1B, FC2A, FC3A, FCRLA, FCRLB) 

#Plot
Plot <- wrap_plots(List,ncol = 3,nrow = 2)

Plot + plot_annotation(tag_levels = 'A')  +
  plot_layout(guides = 'collect')

#Size W1242, H700, TIFF