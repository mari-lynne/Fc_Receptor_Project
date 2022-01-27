library(tidyverse)
#install.packages("tidyverse")#


setwd("C:/Users/mjohnson/Documents/genetics/plink-1.07-dos")
pca <- read_table2("./Plink2.eigenvec", col_names = TRUE)
eigenval <- scan("./Plink2.eigenval")

Sex <- rep(NA, length(pca$IID))
Sex[grep("1", pca$IID)] <- "Male"
Sex[grep("2", pca$IID)] <- "Female"
# location
 Study <- rep(NA, length(pca$IID))
Study[grep("1", pca$IID)] <- "P1"
Study[grep("2", pca$IID)] <- "Tyger"
Study[grep("3", pca$IID)] <- "T1"
Study[grep("4", pca$IID)] <- "T2"

# outcome
Outcome <- rep(NA, length(pca$IID))
Outcome[grep("1", pca$IID)] <- "No Typhoid"
Outcome[grep("2", pca$IID)] <- "Typhoid"

# combine - if you want to plot each in different colours
Sex_Study <- paste0(Sex, "_", Study)

# remake data.frame
pca <- as_tibble(data.frame(pca, Outcome, Study))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = Outcome, shape = Study)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))